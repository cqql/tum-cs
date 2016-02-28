/*
 * Copyright 2013, 2014		Maurice Leclaire <leclaire@in.tum.de>
 *				Stephan M. Guenther <moepi@moepi.net>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation.
 *
 * See COPYING for more details.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <argp.h>
#include <endian.h>
#include <signal.h>
#include <math.h>
#include <inttypes.h>

#include <arpa/inet.h>

#include <net/ethernet.h>

#include <netinet/in.h>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_errno.h>

#include <moep80211/system.h>
#include <moep80211/types.h>

#include <moep80211/modules/moep80211.h>
#include <moep80211/modules/ieee8023.h>

#include "../src/util.h"
#include "../src/list.h"

#include "../src/modules/radio/radiotap.h"

#define EULER_E 2.71828182845904523536028747135266249775724709369995

struct lqe {
	struct list_head list;
	u8 sta[IEEE80211_ALEN];
	u16 last;
	u16 ex_recvd;
	u16 ex_goal;
	bool ex_running;
	struct timespec tslast;
	float a;
	float b;
};

struct list_head lqes_our; // calculated link qualities for links to us
struct list_head lqes_thr; // received link qualities for links from us

struct qinfo {
	u8 sta[IEEE80211_ALEN];
	float a;
	float b;
} __attribute__((packed));

struct qinfos {
	//struct moep_hdr_ext hdr;
	u8 count;
	//struct qinfo i[0];
} __attribute__((packed));


#define MOEP_HDR_BEACON	MOEP_HDR_VENDOR_MIN

struct moep_hdr_beacon {
	struct moep_hdr_ext hdr;
	u16 ex_status;
	struct qinfos qi;
} __attribute__((packed));

static moep_dev_t tap;
static moep_dev_t rad;

static sig_atomic_t _run = 1;

const char *argp_program_version = "ptmbacon 1.0";
const char *argp_program_bug_address = "<michaeli@in.tum.de>";

static char args_doc[] = "IF FREQ";

static char doc[] =
"ptmbeacon - a packet transfer module for moep80211 with beacons\n\n"
"  IF                         Use the radio interface with name IF\n"
"  FREQ                       Use the frequency FREQ [in Hz] for the radio\n"
"                             interface; You can use M for MHz.";

enum fix_args {
	FIX_ARG_IF = 0,
	FIX_ARG_FREQ = 1,
	FIX_ARG_CNT
};

static struct argp_option options[] = {
	{
		.name   = "hwaddr",
		.key    = 'a',
		.arg    = "ADDR",
		.flags  = 0,
		.doc    = "Set the hardware address to ADDR"
	},
	{
		.name   = "ipaddr",
		.key    = 'i',
		.arg    = "ADDR",
		.flags  = 0,
		.doc    = "Set the ip address to ADDR"
	},
	{
		.name   = "mtu",
		.key    = 'm',
		.arg    = "SIZE",
		.flags  = 0,
		.doc    = "Set the mtu to SIZE"
	},
	{
		.name   = "rate",
		.key	= 'r',
		.arg    = "RATE | MCS",
		.flags  = 0,
		.doc	= "Set legacy RATE [r*500kbit/s] or MCS index"
	},
	{
		.name   = "ht",
		.key	= 'h',
		.arg    = "HT",
		.flags  = 0,
		.doc	= "Set HT channel width"
	},
	{
		.name   = "gi",
		.key	= 'g',
		.arg    = "GI",
		.flags  = 0,
		.doc	= "Set GI"
	},
	{
		.name   = "ea",
		.key	= 'e',
		.arg    = "EA",
		.flags  = 0,
		.doc	= "Start full load experiment after e beacons."
	},
	{
		.name   = "el",
		.key	= 'l',
		.arg    = "EL",
		.flags  = 0,
		.doc	= "Number of frames to send at full rate."
	},
	{NULL}
};

static error_t parse_opt(int key, char *arg, struct argp_state *state);

static struct argp argp = {
	options,
	parse_opt,
	args_doc,
	doc
};

static struct config {
	char *rad;
	u8 *hwaddr;
	struct in_addr ip;
	int mtu;
	u64 freq;
	int freq1;
	struct {
		u32 it_present;
		u8 rate;
		struct {
			u8 known;
			u8 flags;
			u8 mcs;
		} mcs;
	} rt;
	u64 moep_chan_width;
	u32 experiment_after;
	u32 experiment_length;
} cfg;

static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
	struct config *config = state->input;
	char *endptr = NULL;
	long long int freq;

	switch (key) {
	case 'a':
		if (!(config->hwaddr = ieee80211_aton(arg)))
			argp_failure(state, 1, errno, 
					"Invalid hardware address");
		break;
	case 'i':
		if (!inet_aton(arg, &config->ip))
			argp_failure(state, 1, errno, "Invalid ip address");
		break;
	case 'm':
		config->mtu = strtol(arg, &endptr, 0);
		if (endptr != NULL && endptr != arg + strlen(arg))
			argp_failure(state, 1, errno, "Invalid mtu: %s", arg);
		if (config->mtu <= 0)
			argp_failure(state, 1, errno,
					"Invalid mtu: %d", config->mtu);
		break;
	case 'r':
		if (config->rt.it_present & BIT(IEEE80211_RADIOTAP_MCS)) {
			config->rt.mcs.known |= IEEE80211_RADIOTAP_MCS_HAVE_MCS;
			config->rt.mcs.mcs = atoi(arg);
		}
		else {
			config->rt.it_present |= BIT(IEEE80211_RADIOTAP_RATE);
			config->rt.rate = atoi(arg);
		}
		break;
	case 'h':
		if (config->rt.it_present & BIT(IEEE80211_RADIOTAP_RATE)) {
			config->rt.it_present &= ~BIT(IEEE80211_RADIOTAP_RATE);
			config->rt.mcs.known |= IEEE80211_RADIOTAP_MCS_HAVE_MCS;
			config->rt.mcs.mcs = config->rt.rate;
			config->rt.rate = 0;
		}
		config->rt.it_present |= BIT(IEEE80211_RADIOTAP_MCS);
		config->rt.mcs.known |= IEEE80211_RADIOTAP_MCS_HAVE_BW;
		if (0 == strncasecmp(arg, "ht20", strlen(arg))) {
			config->rt.mcs.flags |= IEEE80211_RADIOTAP_MCS_BW_20;
			config->moep_chan_width = MOEP80211_CHAN_WIDTH_20;
			break;
		}

		if (strlen(arg) != strlen("ht40*"))
			argp_failure(state, 1, errno, "Invalid HT bandwidth: %s", arg);

		if (0 == strncasecmp(arg, "ht40+", strlen(arg))) {
			config->rt.mcs.flags |= IEEE80211_RADIOTAP_MCS_BW_40;
			config->moep_chan_width = MOEP80211_CHAN_WIDTH_40;
			config->freq1 += 10;
			break;
		}
		else if (0 == strncasecmp(arg, "ht40-", strlen(arg))) {
			config->rt.mcs.flags |= IEEE80211_RADIOTAP_MCS_BW_40;
			config->moep_chan_width = MOEP80211_CHAN_WIDTH_40;
			config->freq1 -= 10;
			break;
		}

		argp_failure(state, 1, errno, "Invalid HT bandwidth: %s", arg);
		break;
	case 'g':
		if (config->rt.it_present & BIT(IEEE80211_RADIOTAP_RATE)) {
			config->rt.it_present &= ~BIT(IEEE80211_RADIOTAP_RATE);
			config->rt.mcs.known |= IEEE80211_RADIOTAP_MCS_HAVE_MCS;
			config->rt.mcs.mcs = config->rt.rate;
			config->rt.rate = 0;
		}
		config->rt.it_present |= BIT(IEEE80211_RADIOTAP_MCS);
		config->rt.mcs.known |= IEEE80211_RADIOTAP_MCS_HAVE_GI;
		if (atoi(arg) == 400)
			config->rt.mcs.flags |= IEEE80211_RADIOTAP_MCS_SGI;
		else if (atoi(arg) != 800)
			argp_failure(state, 1, errno, "Invalid GI: %s", arg);
		break;
	case 'e':
		config->experiment_after = atol(arg);
		break;
	case 'l':
		config->experiment_length = atol(arg);
		break;
	case ARGP_KEY_ARG:
		switch (state->arg_num) {
		case FIX_ARG_IF:
			config->rad = arg;
			break;
		case FIX_ARG_FREQ:
			freq = strtoll(arg, &endptr, 0);
			if (freq < 0)
				argp_failure(state, 1, errno, 
					"Invalid frequency: %lld", freq);
			config->freq = freq;
			config->freq1 += freq;
			break;
		default:
			argp_usage(state);
		}
		break;
	case ARGP_KEY_END:
		if (state->arg_num < FIX_ARG_CNT)
			argp_usage(state);
		break;
	default:
		return ARGP_ERR_UNKNOWN;
	}

	return 0;
}

void sigterm(int sig)
{
	_run = 0;
}

static void prepare_moep80211_radiotap_header(struct moep80211_radiotap *rt)
{
	rt->hdr.it_present = cfg.rt.it_present;
	rt->rate = cfg.rt.rate;
	rt->mcs.known = cfg.rt.mcs.known;
	rt->mcs.flags = cfg.rt.mcs.flags;
	rt->mcs.mcs = cfg.rt.mcs.mcs;

	rt->hdr.it_present |= BIT(IEEE80211_RADIOTAP_TX_FLAGS);
	rt->tx_flags = IEEE80211_RADIOTAP_F_TX_NOACK;
}

static void prepare_moep80211_frame_header(struct moep80211_hdr *hdr)
{
	hdr->frame_control = htole16(IEEE80211_FTYPE_DATA |
			IEEE80211_STYPE_DATA);
	hdr->disc = htole32(MOEP80211_FRAME_DISCRIMINATOR);
	memcpy(hdr->ta, cfg.hwaddr, IEEE80211_ALEN);

	static u16 seq = 0;
	hdr->txseq = seq++;
}

static void taph(moep_dev_t dev, moep_frame_t frame)
{
	struct moep80211_radiotap *radiotap;
	struct moep80211_hdr *hdr;
	struct moep_hdr_pctrl *pctrl;
	struct ether_header ether, *etherptr;
	size_t len;

	if (!(etherptr = moep_frame_ieee8023_hdr(frame))) {
		fprintf(stderr, "ptmbeacon: error: no ether header: %s\n",
							strerror(errno));
		moep_frame_destroy(frame);
		return;
	}
	memcpy(&ether, etherptr, sizeof(ether));

	moep_dev_frame_convert(rad, frame);

	if (!(hdr = moep_frame_moep80211_hdr(frame))) {
		fprintf(stderr, "ptmbeacon: error: no moep80211 header: %s\n",
								strerror(errno));
		moep_frame_destroy(frame);
		return;
	}
	prepare_moep80211_frame_header(hdr);
	memcpy(hdr->ra, ether.ether_dhost, IEEE80211_ALEN);
	memcpy(hdr->ta, ether.ether_shost, IEEE80211_ALEN);

	if (!(pctrl = (struct moep_hdr_pctrl *)moep_frame_add_moep_hdr_ext(frame,
									   MOEP_HDR_PCTRL,
									   sizeof(*pctrl)))) {
		fprintf(stderr, "ptmbeacon: error: cannot add pctrl header: %s\n", strerror(errno));
		moep_frame_destroy(frame);
		return;
	}
	pctrl->type = htole16(be16toh(ether.ether_type));
	if (!moep_frame_get_payload(frame, &len)) {
		fprintf(stderr, "ptmbeacon: error: no payload: %s\n", 
							strerror(errno));
		moep_frame_destroy(frame);
		return;
	}
	pctrl->len = htole16(len);

	if (!(radiotap = moep_frame_radiotap(frame))) {
		fprintf(stderr, "ptmbeacon: error: no radiotap header: %s\n", 
							strerror(errno));
		moep_frame_destroy(frame);
		return;
	}

	prepare_moep80211_radiotap_header(radiotap);

	if (moep_dev_tx(rad, frame)) {
		fprintf(stderr, "ptmbeacon: error: failed to send frame: %s\n", strerror(errno));
		moep_frame_destroy(frame);
		return;
	}

	moep_frame_destroy(frame);
}

static struct lqe * find_lqe(struct list_head* l, u8 *addr) {
	struct list_head *e = 0;
	struct lqe *est = 0;
	list_for_each(e, l) {
		struct lqe *foo = container_of(e, struct lqe, list);
		if(memcmp(foo->sta, addr, IEEE80211_ALEN) == 0)
			est = foo;
	}
	if(!est) {
		est = malloc(sizeof(*est));
		if(!est) {
			fprintf(stderr, "ptmbeacon: malloc fail");
			exit(-1);
		}
		bzero(est, sizeof(*est));
		memcpy(est->sta, addr, IEEE80211_ALEN);
		est->a = est->b = 1;
		list_add(&est->list, l);
	}
	return est;
}
double logp(double a, double b) { return log(a) / log(b); }

double iCDFapproc(double n, double alpha, double beta, double p0) {
  double k = 0.0;
  double S = 0.0;
  double t = 1.0;

  double ab = alpha + beta;
  double an = alpha + n;

  double bound =
      p0 * sqrt(ab / alpha / beta) *
      pow(ab, alpha * (logp(alpha, ab) - 1) + beta * (logp(beta, ab) - 1));
  double upperbound = 4 * n;

  while (S < bound && k < upperbound) {
    double bi = beta + k;
    double anbi = alpha + n + beta + k;
    double P = t * sqrt(anbi / an / bi) *
               pow(anbi, an * (logp(an, anbi) - 1) + bi * (logp(bi, anbi) - 1));
    S += P;
    k += 1;
    t *= (n + k - 1) / (k);
  }

  return k - 1.0;
}

static void adjust_ab(struct lqe * est) {
	double def = 1;
	struct timespec cur, cur2;
	clock_gettime(CLOCK_REALTIME, &cur);
	memcpy(&cur2, &cur, sizeof(cur));
	timespecsub(&cur, &est->tslast);
	double tdiff = cur.tv_sec + (cur.tv_nsec / 1e9f);
	memcpy(&est->tslast, &cur2, sizeof(cur));
	double trust = pow(EULER_E, -tdiff/3);
	//printf(" aB %s: %f->%f, %f->%f\n", ieee80211_ntoa(est->sta), est->a, 1 + trust * (est->a - 1), est->b, 1 + trust * (est->b - 1));
	est->a = def + trust * (est->a - def);
	est->b = def + trust * (est->b - def);
}

static void radh(moep_dev_t dev, moep_frame_t frame)
{
	struct moep80211_hdr *hdr;
	struct moep_hdr_pctrl *pctrl;
	struct moep_hdr_beacon *beacon;
	struct ether_header ether, *etherptr;

	if (!(hdr = moep_frame_moep80211_hdr(frame))) {
		fprintf(stderr, "ptmbeacon: error: no moep80211 header: %s\n", strerror(errno));
		moep_frame_destroy(frame);
		return;
	}

	if (!memcmp(hdr->ta, cfg.hwaddr, IEEE80211_ALEN)) {
		moep_frame_destroy(frame);
		return;
	}

	struct lqe * est = find_lqe(&lqes_our, hdr->ta);
	
	if(est->tslast.tv_sec == 0 && est->tslast.tv_nsec == 0) {
		/* First packet in link */
		clock_gettime(CLOCK_REALTIME, &est->tslast);
		est->last = hdr->txseq;
	} else {
		/* Calculate quality of link packet was received over */
		u64 m = ((u64)hdr->txseq + 0x10000 - est->last) % 0x10000;
		est->last = hdr->txseq;
		adjust_ab(est);
		est->a += 1;
		est->b += m - 1;
	}

	/* Extract received LQs */
	if ((beacon = (struct moep_hdr_beacon *)moep_frame_moep_hdr_ext(frame,
								         MOEP_HDR_BEACON))) {
		if(beacon->qi.count) {
			size_t pll;
			struct qinfo *qis = (struct qinfo*)moep_frame_get_payload(frame, &pll);
			size_t i;
			for(i = 0; i < beacon->qi.count; ++i) {
				struct qinfo *qi = &qis[i];
				if((u8*)qi >= ((u8*)qis) + pll) {
					printf(" warning: evil frame.\n");
					break;
				}
				double pmin = gsl_cdf_beta_Pinv(0.1,qi->a,qi->b);
				u64 l = cfg.experiment_length ? cfg.experiment_length : 64;
				double rmin2 = (iCDFapproc(l, qi->a, qi->b, 0.9) + cfg.experiment_length) / l;
				if(memcmp(cfg.hwaddr, qi->sta, IEEE80211_ALEN) == 0) {
					printf(" qo %s: %f %f -> old: %f new: %f\n", ieee80211_ntoa(qi->sta), qi->a, qi->b, 1/pmin, rmin2);
					struct lqe * est = find_lqe(&lqes_thr, hdr->ta);		
					clock_gettime(CLOCK_REALTIME, &est->tslast);
					est->a = qi->a, est->b = qi->b;
					if(beacon->ex_status != 0 && est->ex_running)
						est->ex_recvd++;
					else if(beacon->ex_status != 0  && !est->ex_running)
						est->ex_running = true,
						est->ex_goal = le16toh(beacon->ex_status),
						est->ex_recvd = 1;
					else if(beacon->ex_status == 0 && est->ex_running)
						printf("Experiment by %s finished: %d / %d\n", ieee80211_ntoa(qi->sta), est->ex_recvd, est->ex_goal), 
						fflush(stdout),
						est->ex_running = false;
				}
			}
		} else
			printf(" no lqi\n");
		moep_frame_destroy(frame);
		return;
	}

	if (!(pctrl = (struct moep_hdr_pctrl *)moep_frame_moep_hdr_ext(frame,
								       MOEP_HDR_PCTRL))) {
		fprintf(stderr, "ptmbeacon: error: no pctrl header: %s\n", strerror(errno));
		moep_frame_destroy(frame);
		return;
	}

	memcpy(ether.ether_dhost, hdr->ra, IEEE80211_ALEN);
	memcpy(ether.ether_shost, hdr->ta, IEEE80211_ALEN);
	ether.ether_type = htobe16(le16toh(pctrl->type));

	if (!moep_frame_adjust_payload_len(frame, le16toh(pctrl->len))) {
		fprintf(stderr, "ptmbeacon: error: failed to adjust payload len: %s\n", strerror(errno));
		moep_frame_destroy(frame);
		return;
	}

	moep_dev_frame_convert(tap, frame);

	if (!(etherptr = moep_frame_ieee8023_hdr(frame))) {
		fprintf(stderr, "ptmbeacon: error: no ether header: %s\n", strerror(errno));
		moep_frame_destroy(frame);
		return;
	}
	memcpy(etherptr, &ether, sizeof(ether));

	if (moep_dev_tx(tap, frame)) {
		fprintf(stderr, "ptmbeacon: error: failed to send frame: %s\n", strerror(errno));
		moep_frame_destroy(frame);
		return;
	}

	moep_frame_destroy(frame);
}

static void send_beacon(u16 ex_status)
{
	moep_frame_t frame;
	struct moep80211_radiotap *radiotap;
	struct moep80211_hdr *hdr;
	struct moep_hdr_beacon *beacon;

	if (!(frame = moep_dev_frame_create(rad))) {
		fprintf(stderr, "ptmbeacon: cannot create frame: %s\n", strerror(errno));
		return;
	}

	if (!(hdr = moep_frame_moep80211_hdr(frame))) {
		fprintf(stderr, "ptmbeacon: error: no moep80211 header: %s\n", strerror(errno));
		moep_frame_destroy(frame);
		return;
	}
	prepare_moep80211_frame_header(hdr);
	memset(hdr->ra, 0xff, IEEE80211_ALEN);
	memcpy(hdr->ta, cfg.hwaddr, IEEE80211_ALEN);

	if (!(beacon = (struct moep_hdr_beacon *)moep_frame_add_moep_hdr_ext(frame,
									     MOEP_HDR_BEACON,
									     sizeof(*beacon)))) {
		fprintf(stderr, "ptmbeacon: error: cannot create beacon header: %s\n", strerror(errno));
		moep_frame_destroy(frame);
		return;
	}

	struct list_head *e = 0;
	size_t s = 0;
	list_for_each(e, &lqes_our) ++s;
	const size_t pll = s*sizeof(struct qinfo);
	struct qinfo *is = malloc(pll); // TODO (actually never do) check alloc.
	size_t i = 0;
	list_for_each(e, &lqes_our) {
		struct lqe *src = container_of(e, struct lqe, list);
		struct qinfo *dst = &is[i];
		memcpy(dst->sta, src->sta, IEEE80211_ALEN);
		adjust_ab(src);
		dst->a = src->a;
		dst->b = src->b;
		++i;
	}
	if(s != 0)
		moep_frame_set_payload(frame, (u8*)is, pll);
	beacon->qi.count = s;
	beacon->ex_status = htole16(ex_status);

	if (!(radiotap = moep_frame_radiotap(frame))) {
		fprintf(stderr, "ptmbeacon: error: no radiotap header: %s\n", strerror(errno));
		moep_frame_destroy(frame);
		free(is);
		return;
	}
	prepare_moep80211_radiotap_header(radiotap);

	if (moep_dev_tx(rad, frame)) {
		fprintf(stderr, "ptmbeacon: error: failed to send frame\n");
	}

	moep_frame_destroy(frame);
	free(is);
}

static int run()
{
	struct timespec interval, tmp, timeout;
	int tx_event;
	sigset_t blockset, oldset;
	fd_set ior;

	if ((tx_event = dup(moep_dev_get_tx_event(rad))) < 0) {
		fprintf(stderr, "ptmbeacon: error: %s\n", strerror(errno));
		free(cfg.hwaddr);
		moep_dev_close(rad);
		moep_dev_close(tap);
		return -1;
	}

	interval.tv_sec = 0;
	interval.tv_nsec = 2e8;

	u32 f_till_experiment = cfg.experiment_after;
	u32 ex_f_count = 0;
	u32 ex_f_max = 0 /* silence a warning */;

	sigfillset(&blockset);
	if (sigprocmask(SIG_SETMASK, &blockset, &oldset)) {
		fprintf(stderr, "ptmbeacon: error: %s\n", strerror(errno));
		free(cfg.hwaddr);
		moep_dev_close(rad);
		moep_dev_close(tap);
		return -1;
	}

	while (_run) {
		FD_ZERO(&ior);
		FD_SET(tx_event, &ior);

		clock_gettime(CLOCK_REALTIME, &timeout);
		if (moep_select(tx_event + 1, &ior, NULL, NULL, NULL, &oldset) < 0) {
			if (errno != EINTR) {
				fprintf(stderr, "ptmbeacon: error: %s\n", strerror(errno));
				free(cfg.hwaddr);
				moep_dev_close(rad);
				moep_dev_close(tap);
				return -1;
			}
		}

		if (!_run)
			break;

		//printf("%d %d %d\n", ex_f_count, f_till_experiment, cfg.experiment_length);
		if((!ex_f_count && f_till_experiment) || !cfg.experiment_length) {
			f_till_experiment--;
			send_beacon(0);

			clock_gettime(CLOCK_REALTIME, &tmp);
			timespecsub(&timeout, &tmp);
			while (timeout.tv_sec < 0)
				timespecadd(&timeout, &interval);
			if (moep_select(0, NULL, NULL, NULL, &timeout, &oldset) < 0) {
				if (errno != EINTR) {
					fprintf(stderr, "ptmbeacon: error: %s\n", strerror(errno));
					free(cfg.hwaddr);
					moep_dev_close(rad);
					moep_dev_close(tap);
					return -1;
				}
			}
		} else {
			if(!ex_f_count) {
				f_till_experiment = cfg.experiment_after;
				// pick a random link (we might want to save which we picked, but oh come on, there's only two boxes in this experiment, right?)
				struct list_head *e = 0;
				size_t s = 0;
				list_for_each(e, &lqes_thr) ++s;
				if(!s) continue;
				s = rand() % s;
				struct lqe *est = 0;
				list_for_each(e, &lqes_thr)
					if(!(s--)) {
						est = container_of(e, struct lqe, list);
						break;
					}	
				adjust_ab(est);
				//double pmin = gsl_cdf_beta_Pinv(0.1,est->a,est->b);
				ex_f_count = ex_f_max = iCDFapproc(cfg.experiment_length, est->a, est->b, 0.9) + cfg.experiment_length;
				//ex_f_count = ex_f_max = 1 / pmin * cfg.experiment_length;
				printf(" ex %s: %d\n", ieee80211_ntoa(est->sta), ex_f_count);
			}
			ex_f_count--;
			send_beacon(ex_f_max);
			//printf("sb: %d\n", ex_f_count);
		}
	}

	sigprocmask(SIG_SETMASK, &oldset, NULL);

	return 0;
}

int main(int argc, char **argv)
{
	gsl_set_error_handler_off ();

	INIT_LIST_HEAD(&lqes_our);
	INIT_LIST_HEAD(&lqes_thr);

	struct sigaction sact;

	sact.sa_handler = sigterm;
	sigfillset(&sact.sa_mask);
	sact.sa_flags = 0;

	if (sigaction(SIGTERM, &sact, NULL)) {
		fprintf(stderr, "ptmbeacon: error: %s\n", strerror(errno));
		return -1;
	}
	if (sigaction(SIGINT, &sact, NULL)) {
		fprintf(stderr, "ptmbeacon: error: %s\n", strerror(errno));
		return -1;
	}

	memset(&cfg, 0, sizeof(cfg));
	cfg.mtu = 1500;
	cfg.moep_chan_width = MOEP80211_CHAN_WIDTH_20_NOHT;
	cfg.experiment_after = 60;
	cfg.experiment_length = 0;
	argp_parse(&argp, argc, argv, 0, 0, &cfg);

	if (!(tap = moep_dev_ieee8023_tap_open(cfg.hwaddr, &cfg.ip, 24,
				       cfg.mtu +
				       sizeof(struct ether_header)))) {
		fprintf(stderr, "ptmbeacon: error: %s\n", strerror(errno));
		return -1;
	}
	if (!(rad = moep_dev_moep80211_open(cfg.rad, cfg.freq,
				    cfg.moep_chan_width,
				    cfg.freq1, 0, cfg.mtu + radiotap_len(-1) +
				    sizeof(struct moep80211_hdr) +
				    sizeof(struct moep_hdr_pctrl)))) {
		fprintf(stderr, "ptmbeacon: error: %s\n", strerror(errno));
		moep_dev_close(tap);
		return -1;
	}

	if (!cfg.hwaddr) {
		if (!(cfg.hwaddr = malloc(IEEE80211_ALEN))) {
			fprintf(stderr, "ptmbeacon: error: failed to allocate memory\n");
			moep_dev_close(rad);
			moep_dev_close(tap);
			return -1;
		}
		if (moep_dev_tap_get_hwaddr(tap, cfg.hwaddr)) {
			fprintf(stderr, "ptmbeacon: error: failed to retrieve hardware address\n");
			free(cfg.hwaddr);
			moep_dev_close(rad);
			moep_dev_close(tap);
			return -1;
		}
	}

	moep_dev_set_rx_handler(tap, taph);
	moep_dev_set_rx_handler(rad, radh);
	moep_dev_pair(tap, rad);

	(void) run();

	free(cfg.hwaddr);
	moep_dev_close(rad);
	moep_dev_close(tap);
	return 0;
}

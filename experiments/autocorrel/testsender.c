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

#include <moep80211/system.h>
#include <moep80211/types.h>

#include <moep80211/modules/moep80211.h>
#include <moep80211/modules/ieee8023.h>

#include "../src/util.h"
#include "../src/list.h"

#include "../src/modules/radio/radiotap.h"

#define MOEP_HDR_BEACON	MOEP_HDR_VENDOR_MIN

struct moep_hdr_beacon {
	struct moep_hdr_ext hdr;
} __attribute__((packed));

struct frinfo {
	u64 seq;
	u64 sec;
	u64 nsec;
} __attribute__((packed));

static u64 seq = 1;

static moep_dev_t rad;
static moep_dev_t tap;

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
		.name   = "br",
		.key	= 'b',
		.arg    = "BR",
		.flags  = 0,
		.doc	= "Set the bacon rate"
	},
	{
		.name   = "al",
		.key	= 'l',
		.arg    = "AL",
		.flags  = 0,
		.doc	= "Pack in some additional dummy payload"
	},
	{
		.name   = "mc",
		.key	= 'c',
		.arg    = "MC",
		.flags  = 0,
		.doc	= "Exit after reaching sequence number 'c'"
	},
	{
		.name   = "bu",
		.key	= 'u',
		.arg    = "BU",
		.flags  = 0,
		.doc	= "Only make the first u frames be beacons, send the rest at maximum rate."
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
	float beaconrate;
	size_t beaconuntil;
	size_t addll;
	size_t maxseq;
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
	case 'b':
		config->beaconrate = atof(arg);
		break;
	case 'l':
		config->addll = atol(arg);
		break;
	case 'c':
		config->maxseq = atol(arg);
		break;
	case 'u':
		config->beaconuntil = atol(arg);
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
	moep_frame_destroy(frame);
}

static void radh(moep_dev_t dev, moep_frame_t frame)
{
	struct moep80211_hdr *hdr;

	if (!(hdr = moep_frame_moep80211_hdr(frame))) {
		fprintf(stderr, "ptmbeacon: error: no moep80211 header: %s\n", strerror(errno));
		moep_frame_destroy(frame);
		return;
	}
	
	if (!memcmp(hdr->ta, cfg.hwaddr, IEEE80211_ALEN)) {
		moep_frame_destroy(frame);
		return;
	}


	size_t pll;
	struct frinfo *is = (struct frinfo*)moep_frame_get_payload(frame, &pll);
	if(pll >= sizeof(struct frinfo)) {
		size_t addl = pll - sizeof(struct frinfo);
		struct timespec tmp;
		clock_gettime(CLOCK_MONOTONIC, &tmp);
		printf("%"PRIu64"-%"PRIu64" %"PRIu64"-%"PRIu64" %zd %"PRIu64"\n", is->sec, is->nsec, tmp.tv_sec, tmp.tv_nsec, addl, is->seq);
	}

	moep_frame_destroy(frame);
}

static void send_beacon()
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

	size_t pll = sizeof(struct frinfo) + cfg.addll;
	struct frinfo *is = malloc(pll);
	struct timespec tmp;
	clock_gettime(CLOCK_MONOTONIC, &tmp);
	is->sec = tmp.tv_sec;
	is->nsec = tmp.tv_nsec;
	is->seq = seq++;
	moep_frame_set_payload(frame, (u8*)is, pll);

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
	//printf("Faia!\n");

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
		return -1;
	}

	interval.tv_sec = 1/cfg.beaconrate;
	interval.tv_nsec = ((int)(1e9/cfg.beaconrate))%1000000000;
	if (cfg.beaconrate > 1e-9)
		fprintf(stderr, "baconing at %"PRIu64"-%"PRIu64".\n", interval.tv_sec, interval.tv_nsec);

	sigfillset(&blockset);
	if (sigprocmask(SIG_SETMASK, &blockset, &oldset)) {
		fprintf(stderr, "ptmbeacon: error: %s\n", strerror(errno));
		free(cfg.hwaddr);
		moep_dev_close(rad);
		return -1;
	}

	while (_run) {
		FD_ZERO(&ior);
		FD_SET(tx_event, &ior);

		if (cfg.beaconrate < 1e9)
			clock_gettime(CLOCK_MONOTONIC, &timeout);
		if (moep_select(tx_event + 1, &ior, NULL, NULL, NULL, &oldset) < 0) {
			if (errno != EINTR) {
				fprintf(stderr, "ptmbeacon: error: %s\n", strerror(errno));
				free(cfg.hwaddr);
				moep_dev_close(rad);
				return -1;
			}
		}

		if (!_run)
			break;

		send_beacon();

		if(seq == cfg.maxseq)
			break;

		if (cfg.beaconrate < 1e9 && ((seq-1) <= cfg.beaconuntil || cfg.beaconuntil == -1)) {
			clock_gettime(CLOCK_MONOTONIC, &tmp);
			timespecsub(&timeout, &tmp);
			while (timeout.tv_sec < 0)
				timespecadd(&timeout, &interval);
			if (moep_select(0, NULL, NULL, NULL, &timeout, &oldset) < 0) {
				if (errno != EINTR) {
					fprintf(stderr, "ptmbeacon: error: %s\n", strerror(errno));
					free(cfg.hwaddr);
					moep_dev_close(rad);
					return -1;
				}
			}
		}
	}

	sigprocmask(SIG_SETMASK, &oldset, NULL);

	return 0;
}

int main(int argc, char **argv)
{
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
	cfg.beaconrate = 20;
	cfg.mtu = 1500;
	cfg.moep_chan_width = MOEP80211_CHAN_WIDTH_20_NOHT;
	cfg.maxseq = -1;
	cfg.beaconuntil = -1;
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

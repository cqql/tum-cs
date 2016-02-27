#!/usr/bin/env python3

import math
from fractions import Fraction

import numpy as np
import scipy.stats as stats


class Timer:
    """Generates a series of timesteps."""

    @staticmethod
    def add_arguments(parser):
        group = parser.add_argument_group("Timers")
        group.add_argument("-t",
                           "--timer",
                           choices=["const", "exp"],
                           default="exp",
                           help="Type of timer to use")

        ConstantStepTimer.add_arguments(parser)
        ExponentialTimer.add_arguments(parser)
        ExperimentTimer.add_arguments(parser)

    @staticmethod
    def from_arguments(arguments):
        if arguments.experiment:
            return ExperimentTimer.from_arguments(arguments)
        elif arguments.timer == "const":
            return ConstantStepTimer.from_arguments(arguments)
        elif arguments.timer == "exp":
            return ExponentialTimer.from_arguments(arguments)

    def generate(self):
        pass


class ConstantStepTimer(Timer):
    """Generates evenly spaced timesteps."""

    @staticmethod
    def add_arguments(parser):
        group = parser.add_argument_group("Constant Step Timer")
        group.add_argument("--step",
                           type=float,
                           default=0.2,
                           help="Length of a timestep")

    @staticmethod
    def from_arguments(arguments):
        return ConstantStepTimer(arguments.time, arguments.step)

    def __init__(self, max, step):
        self.min = 0
        self.max = max
        self.step = step

    def generate(self):
        return np.arange(self.min, self.max, self.step)


class ExponentialTimer(Timer):
    """Draws timesteps from an exponential distribution."""

    @staticmethod
    def add_arguments(parser):
        group = parser.add_argument_group("Exponential Timer")
        group.add_argument("--packet-rate",
                           type=float,
                           default=100,
                           help="Average number of packets per second")

    @staticmethod
    def from_arguments(arguments):
        return ExponentialTimer(arguments.time, arguments.packet_rate)

    def __init__(self, max, rate):
        self.max = max
        self.scale = 1 / rate

    def generate(self):
        t = 0

        while t < self.max:
            dt = np.random.exponential(self.scale)
            t += dt
            yield t


class ExperimentTimer(Timer):
    """Reads timesteps from an experiment log."""

    @staticmethod
    def add_arguments(parser):
        pass

    @staticmethod
    def from_arguments(arguments):
        return ExperimentTimer(arguments.experiment)

    def __init__(self, file):
        self.file = open(file, "r")

        # Drop first line
        self.file.readline()

    def generate(self):
        for line in self.file:
            sec, nano = map(float, line.split(" ")[1].split("-"))

            yield sec + nano * 1e-9


class RateFunction:
    """Computes the true packet success rate at each point in time."""

    @staticmethod
    def add_arguments(parser):
        group = parser.add_argument_group("Rate Functions")
        group.add_argument("-r",
                           "--rate-function",
                           choices=["const", "sin"],
                           default="const",
                           help="Type of rate function to use")

        ConstantRateFunction.add_arguments(parser)
        SineRateFunction.add_arguments(parser)
        UnknownRateFunction.add_arguments(parser)

    @staticmethod
    def from_arguments(arguments):
        if arguments.experiment:
            return UnknownRateFunction.from_arguments(arguments)
        elif arguments.rate_function == "const":
            return ConstantRateFunction.from_arguments(arguments)
        elif arguments.rate_function == "sin":
            return SineRateFunction.from_arguments(arguments)

    def at(self, time):
        pass


class ConstantRateFunction(RateFunction):
    """Takes a constant value everywhere."""

    @staticmethod
    def add_arguments(parser):
        group = parser.add_argument_group("Constant Rate Function")
        group.add_argument("--rate",
                           type=float,
                           default=0.5,
                           help="Packet success rate")

    @staticmethod
    def from_arguments(arguments):
        return ConstantRateFunction(arguments.rate)

    def __init__(self, rate):
        self.rate = max(0, min(1, rate))

    def at(self, time):
        return self.rate


class SineRateFunction(RateFunction):
    """sin(t)."""

    @staticmethod
    def add_arguments(parser):
        group = parser.add_argument_group("Sine Rate Function")
        group.add_argument("--sin-stretch",
                           type=float,
                           default=1.0,
                           help="Stretching factor of the sine curve.")

    @staticmethod
    def from_arguments(arguments):
        return SineRateFunction(arguments.sin_stretch)

    def __init__(self, stretch):
        self.stretch = stretch

    def at(self, time):
        return 0.5 + np.sin(time / self.stretch) / 2


class UnknownRateFunction(RateFunction):
    """The true rate is unknown and therefore NaN."""

    @staticmethod
    def add_arguments(parser):
        pass

    @staticmethod
    def from_arguments(arguments):
        return UnknownRateFunction()

    def at(self, time):
        return np.nan


class Generator:
    """Generates sample packet success and loss."""

    @staticmethod
    def add_arguments(parser):
        group = parser.add_argument_group("Generators")
        group.add_argument("-g",
                           "--generator",
                           choices=["const", "sampling"],
                           default="sampling",
                           help="Type of generator to use")

        ConstantRecvLossGenerator.add_arguments(parser)
        SamplingGenerator.add_arguments(parser)
        ExperimentGenerator.add_arguments(parser)

    @staticmethod
    def from_arguments(arguments):
        if arguments.experiment:
            return ExperimentGenerator.from_arguments(arguments)
        elif arguments.generator == "const":
            return ConstantRecvLossGenerator.from_arguments(arguments)
        elif arguments.generator == "sampling":
            return SamplingGenerator.from_arguments(arguments)

    def sample(self, rate):
        pass


class ConstantRecvLossGenerator(Generator):
    """Receives and loses a constant number of packets at each timestep."""

    @staticmethod
    def add_arguments(parser):
        pass

    @staticmethod
    def from_arguments(arguments):
        return ConstantRecvLossGenerator()

    def sample(self, rate):
        # Approximate rate with a maximum of `total` packets per timestep
        total = 10
        recv = round(total * rate)
        lost = round(total * (1 - rate))

        # Remove common factors
        r = Fraction(recv, lost)

        return (r.numerator, r.denominator)


class SamplingGenerator(Generator):
    """Samples receives and losses according to the given rate."""

    @staticmethod
    def add_arguments(parser):
        pass

    @staticmethod
    def from_arguments(arguments):
        return SamplingGenerator()

    def sample(self, rate):
        if np.random.sample() < rate:
            return (1, 0)
        else:
            return (0, 1)


class ExperimentGenerator(Generator):
    """Read receives and losses from an experiment file."""

    @staticmethod
    def add_arguments(parser):
        pass

    @staticmethod
    def from_arguments(arguments):
        return ExperimentGenerator(arguments.experiment)

    def __init__(self, file):
        self.prevseq = None
        self.file = open(file, "r")

        # Drop first line
        self.file.readline()

    def sample(self, rate):
        line = self.file.readline()
        seq = int(line.split(" ")[-1])

        if self.prevseq is None:
            lost = 0
        else:
            lost = seq - self.prevseq - 1

        self.prevseq = seq

        return (1, lost)


class Estimator:
    """Estimates the packet success rate from delivery and loss data."""

    @staticmethod
    def add_arguments(parser):
        group = parser.add_argument_group("Estimators")
        group.add_argument("-e",
                           "--estimator",
                           choices=["trust", "exp-mean"],
                           default="trust",
                           help="Type of estimator")

        ExpWindowMeanEstimator.add_arguments(parser)
        TrustFunctionEstimator.add_arguments(parser)

    @staticmethod
    def from_arguments(arguments):
        if arguments.estimator == "trust":
            return TrustFunctionEstimator.from_arguments(arguments)
        elif arguments.estimator == "exp-mean":
            return ExpWindowMeanEstimator.from_arguments(arguments)

    def train(self, time, data):
        pass

    def estimate(self, time):
        pass


class ExpWindowMeanEstimator(Estimator):
    """An estimator based on the exponential window mean."""

    @staticmethod
    def add_arguments(parser):
        group = parser.add_argument_group("Exponential Window Mean Estimator")
        group.add_argument("--window-length",
                           type=float,
                           default=0.25,
                           help="Length of the window")
        group.add_argument("--window-initial",
                           type=float,
                           default=0.0,
                           help="Initial guess")
        group.add_argument("--window-rate",
                           type=float,
                           default=0.9,
                           help="Learning rate")

    @staticmethod
    def from_arguments(arguments):
        return ExpWindowMeanEstimator(arguments.window_length,
                                      arguments.window_initial,
                                      arguments.window_rate)

    def __init__(self, length, initial, rate):
        self.length = length
        self.initial = initial
        self.mean = initial
        self.rate = rate
        self.window = None
        self.data = (0, 0)

    def train(self, time, data):
        if self.window is None:
            self.window = time

        self.data = (self.data[0] + data[0], self.data[1] + data[1])

        # Go through as many window as theoretically went by
        while time > self.window + self.length:
            # Update mean and reset
            alpha = self.rate
            if self.data == (0, 0):
                mu = self.initial
            else:
                mu = self.data[0] / (self.data[0] + self.data[1])

            self.mean = alpha * self.mean + (1 - alpha) * mu

            self.window += self.length
            self.data = (0, 0)

    def estimate(self, time):
        return self.mean


class TrustFunctionEstimator(Estimator):
    """An estimator based on the trust function model."""

    @staticmethod
    def add_arguments(parser):
        group = parser.add_argument_group("Trust Function Estimator")
        group.add_argument("--trust-function",
                           choices=["1", "0", "exp"],
                           default="exp",
                           help="Trust function to use")
        group.add_argument("--trust-exp-scale",
                           type=float,
                           default=1.0,
                           help="How fast to forget (greater -> faster)")

    @staticmethod
    def from_arguments(arguments):
        if arguments.trust_function == "1":
            trust = lambda t: 1
        elif arguments.trust_function == "0":

            def trust(t):
                # Trust for a very short time so that you can see something in
                # the plot
                if t < 0.02:
                    return 1
                else:
                    return 0
        elif arguments.trust_function == "exp":
            trust = lambda t: math.exp(-t * arguments.trust_exp_scale)

        return TrustFunctionEstimator(trust)

    def __init__(self, trust):
        self.trust = trust
        self.estimates = []

    def train(self, time, data):
        if len(self.estimates) == 0:
            self.estimates = [(time, 1 + data[0], 1 + data[1])]
        else:
            latest, alpha, beta = self.estimates[-1]
            trust = self.trust(time - latest)
            alpha = 1 + trust * (alpha - 1) + data[0]
            beta = 1 + trust * (beta - 1) + data[1]

            self.estimates.append((time, alpha, beta))

    def estimate(self, time):
        if len(self.estimates) == 0:
            return stats.beta(1, 1)
        else:
            latest, alpha, beta = self.estimates[-1]
            trust = self.trust(time - latest)
            alpha = 1 + trust * (alpha - 1)
            beta = 1 + trust * (beta - 1)

            return stats.beta(alpha, beta)


class Simulation:
    @staticmethod
    def add_arguments(parser):
        group = parser.add_argument_group("Simulation")
        group.add_argument("--sampling-rate",
                           type=float,
                           default=100,
                           help="How often to sample per second")
        group.add_argument("--start",
                           type=float,
                           default=0,
                           help="Start of plot in normalized time")
        group.add_argument("--end",
                           type=float,
                           default=None,
                           help="End of plot in normalized time")

    @staticmethod
    def from_arguments(arguments):
        return Simulation(arguments.sampling_rate, arguments.start,
                          arguments.end)

    def __init__(self, sampling_rate, start, end):
        self.sampling_rate = sampling_rate
        self.start = start
        self.end = end

    def run(self, timer, rate, generator, estimators):
        measurements = []
        estimations = [(e, []) for e in estimators]
        minT = None
        t = 0
        timestep = 1.0 / self.sampling_rate

        for time in timer.generate():
            if minT is None:
                minT = time

            # Normalize times to start from 0
            time -= minT

            finaliteration = False
            if self.end and time > self.end:
                time = self.end
                finaliteration = True

            # Sample estimates until the next data point
            while t < time:
                if self.start <= t:
                    for e, data in estimations:
                        data.append((t, rate.at(t), e.estimate(t)))

                t = t + timestep

            # Train with next data point
            sample = generator.sample(rate.at(time))

            if self.start <= t:
                measurements.append((time, *sample))

            for e, data in estimations:
                e.train(time, sample)

            if finaliteration:
                break

        return (measurements, estimations)

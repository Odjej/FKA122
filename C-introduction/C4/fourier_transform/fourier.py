
import numpy as np
import matplotlib.pyplot as plt
import argparse


# --------------------------------------------------
# Signal generation
# --------------------------------------------------

def generate_signal(a=1.0, f=2.0, phi=0.0, dt=0.1, N=250):
    """
    Generate h(t) = a*cos(2*pi*f*t + phi).

    Returns:
        t: Time array
        h: Signal array
    """
    t = np.arange(N) * dt
    h = a * np.cos(2 * np.pi * f * t + phi)

    return t, h


def generate_two_signals(a1=1.0, a2=1.0,
                         f1=2.0, f2=6.0,
                         phi1=0.0, phi2=0.0,
                         dt=0.1, N=250):
    """
    Generate the sum of two cosine signals.
    """
    t = np.arange(N) * dt

    h = (a1 * np.cos(2 * np.pi * f1 * t + phi1)
         + a2 * np.cos(2 * np.pi * f2 * t + phi2))

    return t, h


# --------------------------------------------------
# Fourier transform and power spectrum
# --------------------------------------------------

def calculate_spectrum(h, dt):
    """
    Calculate the DFT and asymmetric power spectrum.

    Returns:
        frequencies_asymmetric: Asymmetric frequency array
        P: Power spectrum
        frequencies_shifted: Shifted frequency array
        P_shifted: Shifted power spectrum
    """
    N = len(h)

    H = np.fft.fft(h)
    P = np.abs(H)**2 / N
    
    # Asymmetric frequency axis: 0 to approximately fs
    frequencies_asymmetric = np.arange(N) / (N * dt)

    # Symmetric frequency axis: negative to positive
    frequencies = np.fft.fftfreq(N, d=dt)
    frequencies_shifted = np.fft.fftshift(frequencies)
    P_shifted = np.fft.fftshift(P)

    return frequencies_asymmetric, P, frequencies_shifted, P_shifted


# --------------------------------------------------
# Plotting
# --------------------------------------------------

def plot_signal(t, h, title="Signal"):
    plt.figure(figsize=(9, 4))
    plt.plot(t, h)
    plt.xlabel("Time")
    plt.ylabel("h(t)")
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def plot_spectrum(frequencies, P, title="Power spectrum",
                  symmetric=False):
    plt.figure(figsize=(9, 4))

    plt.plot(frequencies, P)
    plt.xlabel("Frequency")
    plt.ylabel("Power spectrum")
    plt.title(title)

    plt.grid(True)
    plt.tight_layout()
    plt.show()


# --------------------------------------------------
# Experiments
# --------------------------------------------------

def experiment_signal(args):
    """Plot a single cosine signal."""
    t, h = generate_signal(
        a=args.amplitude,
        f=args.frequency,
        phi=args.phase,
        dt=args.dt,
        N=args.N
    )

    plot_signal(t, h)


def experiment_spectrum(args):
    """Plot both asymmetric and symmetric power spectra."""
    t, h = generate_signal(
        a=args.amplitude,
        f=args.frequency,
        phi=args.phase,
        dt=args.dt,
        N=args.N
    )

    frequencies, P, f_shifted, P_shifted = (
    calculate_spectrum(h, args.dt)
    )

    if args.signal_plot:
        plot_signal(t, h)

    if args.symmetric:
        plot_spectrum(
            f_shifted,
            P_shifted,
            title="Symmetric power spectrum"
        )
    else:
        plot_spectrum(
            frequencies,
            P,
            title="Asymmetric power spectrum"
        )


def experiment_two_frequencies(args):
    """Investigate two frequencies and the Nyquist limit."""
    t, h = generate_two_signals(
        a1=1.0,
        a2=1.0,
        f1=args.f1,
        f2=args.f2,
        dt=args.dt,
        N=args.N
    )

    frequencies, P, f_shifted, P_shifted = (
        calculate_spectrum(h, args.dt)
    )

    nyquist = 1 / (2 * args.dt)

    print(f"Sampling frequency: {1/args.dt:.3f}")
    print(f"Nyquist frequency:  {nyquist:.3f}")
    print(f"Minimum frequency: {-nyquist:.3f}")
    print(f"Maximum frequency: {nyquist - 1/(args.N*args.dt):.3f}")

    plot_signal(t, h, title="Two-frequency signal")

    plot_spectrum(
        f_shifted,
        P_shifted,
        title="Two-frequency symmetric power spectrum"
    )



def main():
    print("Starting Fourier program...")
    
    parser = argparse.ArgumentParser(
        description="Discrete Fourier transform experiments"
    )

    parser.add_argument(
        "--mode",
        choices=["signal", "spectrum", "two"],
        default="spectrum",
        help="Experiment to run"
    )

    parser.add_argument("--N", type=int, default=250)
    parser.add_argument("--dt", type=float, default=0.1)

    parser.add_argument("--amplitude", type=float, default=1.0)
    parser.add_argument("--frequency", type=float, default=2.0)
    parser.add_argument("--phase", type=float, default=0.0)

    parser.add_argument("--f1", type=float, default=2.0)
    parser.add_argument("--f2", type=float, default=6.0)

    parser.add_argument(
        "--symmetric",
        action="store_true",
        help="Plot the symmetric power spectrum"
    )

    parser.add_argument(
        "--signal-plot",
        action="store_true",
        help="Also plot the time-domain signal"
    )

    args = parser.parse_args()

    if args.mode == "signal":
        experiment_signal(args)

    elif args.mode == "spectrum":
        experiment_spectrum(args)

    elif args.mode == "two":
        experiment_two_frequencies(args)


if __name__ == "__main__":
    main()

# The original FFT stores the frequency components in an order where zero frequency comes first, 
# followed by positive frequencies and then negative frequencies.

# The `np.fft.fftshift` function is used to shift the zero frequency component to the center of the spectrum,
# which is often more convenient for visualization and analysis. This is especially useful when plotting the power spectrum, 
# as it allows for a more intuitive representation of the frequency components, with negative frequencies on the left and positive frequencies on the right.

# N = 253 sompletes 50.6 periods, compared to 50 and 51 periods for N = 250 and N = 255, respectively. 
# The signal does not complete an integer number of periods, which can lead to spectral leakage in the Fourier transform.
# We see this in the PS as broadening of the peak around 2 Hz.

# For the two-frequency signal, we first have f_s = 1/dt = 10 Hz, which means the Nyquist frequency is f_N = f_s/2 = 5 Hz.
# The 6 Hz component exceeds the Nyquist frequency, which will lead to aliasing. 
# The 6 Hz component will be reflected back into the frequency range below the Nyquist frequency, 
# resulting in a distorted representation of the original signal in the frequency domain.
# We instead see a peak at 4 Hz, which is the aliased frequency of the 6 Hz component.
# Improve sampling by increasing the sampkling frequency, dt = 0.05, means f_s = 20 Hz, f_N = 10 Hz, 
# and both frequency components are now below the Nyquist frequency, allowing for accurate representation in the frequency domain.
# We then get the correct peaks at 2 Hz and 6 Hz in the power spectrum, accurately representing the original signal's frequency content.

# Nyquist sampling criterion: f_s >= 2 f_max, where f_max is the highest frequency component in the signal.
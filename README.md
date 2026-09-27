# Vital Signs Sensing for 6G Networks

Integrated Sensing and Communication (ISAC) simulation in MATLAB: an OFDM communication link and an FMCW radar sensing link share the same time-frequency resource to simultaneously transmit data and detect a human target's breathing and heartbeat.

## System Overview

- **Communication link**: QPSK-modulated OFDM (configurable subcarriers/CP), with convolutional coding (rate 1/2, `poly2trellis(7,[171 133])`) + hard-decision Viterbi decoding on the receive side.
- **Sensing link**: FMCW radar at 77 GHz carrier / ~4 GHz bandwidth, simulating a target at 2 m range whose round-trip delay is phase-modulated by chest-wall motion from breathing (~0.333 Hz) and heartbeat (~1.2 Hz).
- **Time-division multiplexing**: communication and sensing waveforms share the channel and are separated on the receive side before independent processing.

## Vital Sign Extraction Pipeline

1. Range-FFT on the FMCW receive matrix to locate the target range bin.
2. Phase extraction + unwrapping at that range bin across slow time (chirps).
3. Band-pass filtering to separate the breathing band (0.1–0.5 Hz) from the heartbeat band (0.8–2 Hz).
4. Frequency estimation on each band via three independent methods, compared side by side:
   - **2D-FFT** (`main.m`) — baseline spectral peak picking.
   - **MUSIC** (`MUSIC.m`) — subspace-based super-resolution frequency estimation.
   - **VMD** (`VMD.m`, Variational Mode Decomposition) — adaptively decomposes the signal into modes before estimating each mode's center frequency; best resolution of the three on closely-spaced/weak components.

## Files

| File | Purpose |
|---|---|
| `main.m` | Full ISAC pipeline: OFDM tx/rx + FMCW tx/rx + BER + breathing/heartbeat estimation via 2D-FFT, MUSIC, and VMD. Entry point. |
| `comm.m` | Standalone communication-link Monte Carlo sweep: BER vs. SNR, with and without convolutional coding. |
| `fmcw_breathe_and_heartbeat.m` | Minimal standalone FMCW breathing/heartbeat simulation, isolated from the OFDM link. |
| `simulation_of_breathe_and_heartbeat.m` | Toy signal (two superimposed sinusoids at the target breathing/heartbeat frequencies) used to sanity-check the FFT-based estimation before running it on the full radar signal. |
| `MUSIC.m` | MUSIC pseudo-spectrum frequency estimator, used for vital-sign frequency estimation. |
| `VMD.m` | Variational Mode Decomposition, used to separate and estimate breathing/heartbeat components. |

## Requirements

MATLAB with the Communications Toolbox and Signal Processing Toolbox (for `pskmod`/`pskdemod`, `poly2trellis`/`convenc`/`vitdec`, `awgn`, `butter`/`filter`, `pwelch`).

## Running

Open `main.m` in MATLAB and run it — it prints the coded-link BER and the breathing/heartbeat frequency estimates from all three methods to the console, and plots the range-Doppler map, phase traces, and per-method spectra.

Run `comm.m` separately to reproduce the BER-vs-SNR curve for the communication link alone.

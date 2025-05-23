# Adaptive Equalizer & Notch Filter System

This project implements key adaptive signal processing algorithms in MATLAB to suppress interference and recover distorted signals in communication channels.

## 🔧 Algorithms
- **LMS Equalizer**: Trained with known input-output pairs to adapt to channel distortion
- **CMA Equalizer**: Blind equalization without training sequence
- **Adaptive Notch Filter**: Tracks and removes a dynamic sinusoidal interferer

## 📈 Results
- Verified convergence under varying filter orders and SNRs
- Frequency-domain plots for notch depth and equalizer stability
- Works on 100k-sample data streams

This project was part of ESE 5310 (Digital Signal Processing) at the University of Pennsylvania.

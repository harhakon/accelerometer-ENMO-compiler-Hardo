Hardo compiler is a C based tool used to process accelerometer raw text data files (.csv) with the aim of translating them into an epoch file (.csv) that contains desired time series such as ENMOs, wrist angle, length of nonwear time, MADs and steps.

Works well with the Axivity model accelerometer, other models may need adjustment, especially with calibration coefficients.

- Window peak detection algorithm used to count steps. (Femiano, R. 2022).

- Wrist angle was estimated with formula tan-1(az/sqrt((ax2+ay2))*180*pi   (Van Hees 2015).

- Possible to use Butterworth filtering.

- the length of the nonwear time is determined from the first time point of the period.

- Epoch file format (time, enmo (mad), wrist angle, nonwear time, steps)

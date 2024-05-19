# Source Factor Extraction

My research on extracting the excitation and source factors of binary black hole mergers.

## Excitation Coefficient

The Excitation Coefficient is obtained from a fitting function in the [SourceFactors notebook](https://github.com/Aditya7s/SourceFactors/blob/master/SourceFactors.ipynb) up to the seventh overtone using the quasinormal modes. The values are saved as a .csv file.
The goal of this code is to import SXS data, limiting to spin $\chi = 0.69 \pm 0.05$ or close to it, then fit the $l=2, m=2$ waveform to a damped sinusoid model. After fitting, the coefficients for up to overtone $N=7$ is extracting and stored.

### The Fitting

A non-linear least-squares minimization is done on the waveform and fitted to $$h_{lm}^N(t) = d + \sum_{n=0}^N C_{lmn}e^{-i\omega_{lmn}(t-t_0)}$$ where $h_{lm}^N(t)$ will be the complex strain and $t$ is the time from the waveform in the NR simulation. The coefficients $C_{lmn}$ is what is going to be extracted from the fitting.

The function takes in the start time $T$ which tells how far after the peak strain time should the fitting start $t_0$, then it takes the waveform for the y values in the fitting, then the spin for calculating the qnm frequencies which after dividing after the next input, the mass, will provide the $\omega_{lmn}$. Once the fitting is done until $t_0 + 90M$ the lmfit parameters are returned.

## Excitation and Source Factors

The excitation and source factors are obtained within the [MathematicaNB](https://github.com/Aditya7s/SourceFactors/tree/master/MathematicaNB), specifically the [ExcitationFactors notebook](https://github.com/Aditya7s/SourceFactors/blob/master/MathematicaNB/ExcitationFactors.nb) which utilizes the KerrMST packages, which is a solver for the Teukolsky equations and uses the Mano-Suzuki-Takasuki method to obtain the Excitation Factors. The Source Factors are obtained by dividing the Excitation Coefficients by the Excitation Factors.

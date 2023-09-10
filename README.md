# Source Factor Extraction

My research on extracting the excitation and source factors of binary black hole mergers.

## Excitation Coefficient

The Excitation Coefficient is obtained from a fitting function in the [SourceFactors notebook](SourceFactors.ipynb) up to the seventh overtone using the quasinormal modes. The values are saved as a .csv file.

## Excitation and Source Factors

The excitation and source factors are obtained within the [MathematicaNB](MathematicaNB), specifically the [ExcitationFactors notebook](MathematicaNB/ExcitationFactors.nb) which utilizes the KerrMST packages, which is a solver for Leavers equations and uses the Mano-Suzuki-Takasuki method to obtain the Excitation Factors. The Source Factors are obtained by dividing the Excitation Coefficients by the Excitation Factors.

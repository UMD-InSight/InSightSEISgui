# InSight_seisGUI
A graphical user interface to download, process, and visualise InSight (VBB) seismic data

## Dependencies

The script has some dependencies and the following packages should be installed in order to run it: 
  csv, matplotlib, numpy, obspy, os pillow, requests, scipy, shutil, xml

## Contents

This repository contains:

1. The 12th version of the Marsquakes Catalog. Citation: InSight Marsquake Service (2022). Use the link for the download and cite accordingly.
2. The ELYSE station dataless file. Citation: InSight MARS SEIS Data Service (2019). Use the link for the download and cite accordingly.
3. The InSight_datagui.py script.

## User guidelines

Once all the package dependencies are satisfied, there is no need to install any additional package, just run the script with a python3 command:

```
python3 InSight_datagui.py
```

The graphical interface will open. Then, the following options are available:

- On the left part of the window, there is a selection tool for event quality and frequency category. By selecting filtering the quality type and frequency family the filtered seismic catalog is generated in the window.
- Select the event to download and process and push the button at the bottom of the event list. The data will be downloaded in a "DATA" directory, within your working directory in .mseed format. Raw, displacement, velocity and acceleration data will be downloaded.

Once the "Download and process" button is not anymore pushed, the data are downloaded and you can work with their visualisation.

In order to select the correct frequency band for the seismogram, first generate the spectrogram for the raw, displacement, or acceleration data, by selecting the respective option button and pushing the button "Make spectrograms".

Once the spectrograms are generated, you can save the figure by pushing the "Export figure" button. This button will export the exact figure that you see at the moment that you push the button in a directory named "Figures/Spectrograms".

In order to visualise the seismograms you should define the low and high frequency for the filtering in the first two boxes of the bottom part of the window, as well as the minimum and maximum time, given in seconds after the origin time. Then you push the "Make seismograms" button to visualise the seismograms. As long as the "amplitude" selection is set to 0, the amplitude in the figure is taken automatically from the maximum amplitude of the waveforms. In order to set your own and eventually avoid the representation of glitches amplitude, write the absolut value of the minimum and maximum amplitude that you want to visualise, then push the "Make seismograms" button again.

When you make new seismograms, the spectrograms are showing the minimum and maximum frequency, as well as the minimum and maximum time, with white dashed lines. You can use this feature to select correctly the event associated signal.

When you have the desired figures, you can export both spectrograms and seismograms with the respective "Export figure buttons". The figures are saved in the "Figures/Spectrograms" and "Figures/Seismograms" directories respectively.

### Screenshot

![Alt text](/InSight_datagui_screenshot.png "InSight seisGUI screenshot")


## Citations and acknowledgements

The scripts of this repository download Mars InSight seismic data from IRIS. It is developped by Foivos Karakostas, Ross Maguire, Doeyeon Kim, Nick Schmerr, Ved Lekić, Jessica Irving. Please acknowledge.

When you use InSight SEIS Data, please follow the [citation instructions](https://www.seis-insight.eu/en/science/seis-data/seis-citation-information) that are also copied here:

*SEIS data must be cited as reference in the following way:*

*Citation in text :*

*InSight Mars SEIS Data Service. (2019). SEIS raw data, Insight Mission. IPGP, JPL, CNES, ETHZ, ICL, MPS, ISAE-Supaero, LPG, MFSC. https://doi.org/10.18715/SEIS.INSIGHT.XB_2016*

*In addition an acknowledgement must also be provided to the SEIS operators as follows:*

*"We acknowledge NASA, CNES, their partner agencies and Institutions (UKSA, SSO, DLR, JPL, IPGP-CNRS, ETHZ, IC, MPS-MPG) and the flight operations team at JPL, SISMOC, MSDS, IRIS-DMC and PDS for providing SEED SEIS data."*

*Furthermore, the SEIS experiment paper (Lognonné et al., 2019) must be used as reference for describing the instrument, in addition to the SEIS team papers used by the user for the analysis. For a collection of SEIS papers, see Banerdt and Russel, 2017 and Banerdt and Russel, 2019.*

## References

1. Clinton, J. F., Ceylan, S., van Driel, M., Giardini, D., Stähler, S. C., Böse, M., … Stott, A. E. (2021). The Marsquake catalogue from InSight, sols 0–478. Physics of the Earth and Planetary Interiors, 310, 106595. https://doi:10.1016/j.pepi.2020.106595
2. InSight Mars SEIS Data Service. (2019). SEIS raw data, Insight Mission. IPGP, JPL, CNES, ETHZ, ICL, MPS, ISAE-Supaero, LPG, MFSC. https://doi.org/10.18715/SEIS.INSIGHT.XB_2016
3. InSight Marsquake Service (2022). Mars Seismic Catalogue, InSight Mission; V12 2022-10-01. ETHZ, IPGP, JPL, ICL, MPS, Univ. Bristol. https://doi.org/10.12686/a18
4. Lognonné, P., Banerdt, W.B., Giardini, D. et al. (2019). SEIS: Insight’s Seismic Experiment for Internal Structure of Mars. Space Sci Rev 215, 12. https://doi.org/10.1007/s11214-018-0574-6

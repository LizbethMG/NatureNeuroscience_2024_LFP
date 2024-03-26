<!-- Section 1: Logo and short code description-->  <!--  -->
<a name="readme-top"></a>
<br /> <!-- break -->
<div align="center">  <!-- A block of content centred for the logo and description -->
  <a href="https://github.com/LizbethMG/NatureNeuroscience_2024_LFP">  <!-- CLickable logo; TODO: change to the article link -->
    <img src="images/CL_logo_LFP.png" alt="Logo" width="80" height="100">
  </a>

  <h2 align="center">Mondragon-Gonzalez et al., 2024,  Nature Neuroscience</h2>  <!-- Header tag -->
  <h3 align="center">LFP analysis</h3>  <!-- Header tag -->
  
  <p align="center">
    This repository contains basic local field potential analysis scripts and functions. 
    <br />
    <!-- TODO: change to the article link -->
    <a href="https://github.com/LizbethMG/2024_Mondragon-Gonzalez_NatureNeuroscience"><strong> Link to article (TODO) »</strong></a>
    <br />
  </p>
</div>

## Table of Contents 

- [1 System Requirements](#1-system-requirements)
- [2 Installation Guide](#2-installation-guide)
- [3 Description of Code's Functionality](#3-description-of-codes-functionality)
- [4 Instructions of Use](#4-instructions-of-use)
- [5 Citation](#5-citation)
- [6 License](#6-license)

## 1 System Requirements
* <b> MATLAB version: </b> The code was developped using Matlab R2017b and also tested using Matlab R2019b and R2023b.
* <b> Operating System: </b> The code was developped and tested under Windows 7 and 10.
* <b> MATLAB Toolboxes: </b> Curve Fitting Toolbox, Signal processing toolbox, Wavelet Toolbox ( and helperCWTTimeFreqPlot.m from R2016a, included in the repository)
* <b> Third-party MATLAB functions: </b> Kelly Kearney (2024). ['boundedline.m'](https://github.com/kakearney/boundedline-pkg), GitHub. Retrieved March 25, 2024. Brandon Kuczenski (2024). ['vline.m'](https://www.mathworks.com/matlabcentral/fileexchange/1039-hline-and-vline), MATLAB Central File Exchange. Retrieved March 25, 2024.

<div align="right">[ <a href="#readme-top">↑ Back to top ↑</a> ]</div>

## 2 Installation Guide
<ol>

  <li> <a href="https://docs.github.com/fr/repositories/creating-and-managing-repositories/cloning-a-repository"> <strong> Clone </strong> </a> the repository to your local machine:
    <pre><code>git clone https://github.com/LizbethMG/NatureNeuroscience_2024_LFP.git</code></pre>
  </li>
  <li>Navigate to the cloned directory:
    <pre><code>cd NatureNeuroscience_2024_LFP</code></pre>
  </li>
  <li> Download the data used for the analysis and add it to your 'Data' project's folder. </li>
  <li> Add the project folder to your MATLAB path. You can do this in two ways:</li>

<h5>Using the MATLAB Command Window</h5>
    <ul>
        <li>Open MATLAB.</li>
        <li>In the Command Window, type the following command, replacing <code>yourpath</code> with the path to your project:
            <pre><code>addpath('yourpath');
savepath;</code></pre> </li>
    </ul>
<h5>Using the MATLAB Set Path Dialog</h5>
    <ul>
        <li>Open MATLAB.</li>
        <li>On the Home tab, in the Environment section, click Set Path.</li>
        <li>Click Add with Subfolders.</li>
        <li>Browse to the project folder and click OK.</li>
        <li>Click Save and then Close to save the changes.</li>
    </ul>
    
  <li> Install (if not already done) the Matlab Toolboxes mentioned in <b> "1 System Requirements" </b> </li>
</ol>
Typical install time:  10 min.
<div align="right">[ <a href="#readme-top">↑ Back to top ↑</a> ]</div>

## 3 Description of Code's Functionality

<br> &rarr;  Basic analysis scripts to explore LFP signals around grooming events. The complete description can be found in the Methods section of the paper.
<br> &diams; For more detailed information, users can refer to the in-code comments within the script and its functions.

### INPUT 
*  `config.m`  configuration file of the experiment(s)
*  data files copied in the [Data/](data/) folder 
### Data processing Wavelet analysis 
* `slmg_waveletn.m`
  * Loads data for each experiment, and prepares it for further analysis. It ensures that the data is standardized across events and sessions to facilitate consistent analysis.
  * Low pass filtering
  * Wavelet analysis
  * Compute power within a specified frequency band
  * Compute maximum wavelet power within specified frequency band
#### Visualization
* `plotWaveletContour`  Creates a contour plot of continuous wavelet transform (cwt) coefficients over time and frequency. 
### OUTPUT
* Generates different plots within a frequency band:
  * Wavelet contour 
  * Power vs time
  * Power vs frequency
* `maxPower` Maximum power and its corresponding time and frequency values

<div align="right">[ <a href="#readme-top">↑ Back to top ↑</a> ]</div>

## 4 Instructions for Use

#### Working with the demo data
- Ensure you added the project folder to your MATLAB path and set it as the current folder. In Matlab Command Window type `pwd` the result should be `'YOUR_PATH\NatureNeuroscience_2024_LFP'`
- Run the Live Script file `LFP_Analysis.m`. It will run the code for the data folder selected in `config.m` line  `confWV.animalList = {M#ID};` .
#### Working with Dataset 2
 It features the electrophysiological recordings in the mice's orbitofrontal cortex around grooming events. 
The data is available at the [Open Science Framework](https://osf.io/kdmjt/) **DOI: 10.17605/OSF.IO/KDMJT**

- [Download](https://osf.io/kdmjt/) dataset 2.
- Add the downloaded folder to your MATLAB path containing the code inside the `Data\extractedData` folder.
- To analyze data from subject 1 you would have the source data file as:  `YOUR_PATH\NatureNeuroscience_2024_LFP\Data\extractedData\M1\datasource_M1.mat`   
Each sub-folder in the data set 1 corresponds to a single  subject (SAPAP3-KO). Each experiment has its configuration information in the `config.m`

<b> &diams; Expected run time for demo: </b> For one subject  (20 grooming events, 4 channels) like the one on dataset2 (M1), on a general-purpose computer expect the following approximately runtimes:
<br> ~ 14 s -- CPU: Intel Core i9-7960X @ 2.8 GHz, RAM: 32 GB, Operating System: Windows 10 64-bit, Matlab 2019b.

<div align="right">[ <a href="#readme-top">↑ Back to top ↑</a> ]</div>

## 5 Citation
If you use this code or data we kindly ask you to cite our work. 

- <b> Data: </b>
> (APA style) Mondragón-González, S. L. (2024, March 26). 2024_Mondragon-Gonzalez_NatureNeuroscience. https://doi.org/10.17605/OSF.IO/KDMJT

- <b> Article: </b> Mondragon et al 2024: [TODO: DOI here](https://github.com/LizbethMG/Mondragon_NatNeuro_CL)
> @article{Mondragon2024,
        title = {Closed-loop recruitment of striatal interneurons prevents compulsive-like grooming behaviours},
        author = {Sirenia Lizbeth Mondragón-González and Christiane Schreiweis and Eric Burguière},
        journal = {Nature Neuroscience},
        year = {2024},
        url = { TODO }}

<div align="right">[ <a href="#readme-top">↑ Back to top ↑</a> ]</div>

## 6 License
Copyright <2024> <COPYRIGHT Sirenia Lizbeth Mondragón-González>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
<div align="right">[ <a href="#readme-top">↑ Back to top ↑</a> ]</div>

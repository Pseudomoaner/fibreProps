# fibreProps
Automated detection and measurement of peptidoglycan fibres in AFM images

## What you will need to do before beginning analysis:
1. **Exporting AFM data:**
2. **File structure:** fibreProps assumes that your AFM data is collected into directories that represent different conditions. 'Conditions' here refers to any set of experimental parameters that may be interesting to compare between such as different species identities, different drug treatments, different cell wall regions etc. Specifically, your data should be arranged as:
```
   Root/--Condition1/--Image1.txt
     |         |-------Image2.txt
     |         |-------Image3.txt
     |         :
     |
     |----Condition2/--Image1.txt
     |         |-------Image2.txt
     |         |-------Image3.txt
     |         :
     :
```
Note that '/Condition1/', 'Image1' etc. can be replaced with arbitrary, informative labels.

3. **Defining cell orientations:**
It is often useful to measure fibre orientations relative to some external angle, such as the long axis of the cell. fibrePlotter supports alignment of fibre orientation measurements to such external measurements, but must be told what they are. The CellOrentations.csv file is the means by which this data is provided to fibrePlotter. This file simply consists of Name,Value pairs which indicate the name of each image and the associated angle (in degrees, measured using standard polar coordinate conventions), respectively. A separate CellOrientations.csv file must be generated for each condition's directory.

Examples of the required AFM data format, file structure and CellOrientations.csv files are provided in the \example\ directory.

4. **Installation:** Both GUIs are implemented as Matlab Apps. To install, simply double-click on the relevant .mlappinstall file with Matlab open and follow the on-screen instructions. The app will then appear in the 'Apps' tab:
<p align="center">
  <img src="https://raw.githubusercontent.com/Pseudomoaner/fibreProps/main/Images/ML_Apps.png" alt="Matlab Apps"/>
</p>

## Detecting fibres: FibreFinder

FibreFinder is a Graphical User Interface (GUI) that 

## Plotting fibres: FibrePlotter

FibrePlotter is a separate GUI that enables rapid statistical and visual comparisons between multiple conditions. 

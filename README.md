**Please read carefully the following 3 sections in order to successfully run the model-driven motion correction algorithm.**

## Intallation and Requirements

You will need to have Python (>=3.5) installed in your machine, regardless of its operative system because it's cross-platform ([Download](https://www.python.org/downloads/)).

Then, open a `Terminal` (MacOS and Linux) or a `cmd` Windows and run the following command to install the required packages:

`pip install -r requirements.txt`

Finally, you must have `elastix` installed in your machine in order to perform the motion correction steps in the process. After installation, you must add the downloaded executable to the `$PATH` environment variable of your user/machine, otherwise Python will not detect `elastix`.

Links:

[Elastix Main Page](https://elastix.lumc.nl/)
[Elastix Download](https://github.com/SuperElastix/elastix/releases/tag/5.0.1)


## Schematic of folders structure

`GroupWise_VS_MDR.py` -> main python script to run in order to perform the comparative study and to generate the results described in the paper.

`Github Repository` -> parent folder
                  Inside it are stored all the patient data that will be processed, each patient in its own folder.
                  Inside it are also stored the folders: Elastix_Parameters_Files and AIFs as explained by the schematic below.
                  The path for this folder is given in line 1235 of "GroupWise_VS_MDR.py" as value for the global variable DATA_PATH and this can be changed

.
+-- AIFs (this folder is needed for DCE data and contains the aif for each patient in a .txt form)
|    |
|    +-- patient_001
|    |          |
|    |          +-- AIF__2CFiltration__Curve.txt
|    |
|    +-- patient_002
|    |           |
|    |           +-- AIF__2CFiltration__Curve.txt
|    |  
|    +-- ...etc...
|
|
|
|
+-- Elastix_Parameters_Files (the full path for this folder is given in line 1297 in the argument parameters_file=...)
|    |
|    +-- BSplines_DCE.txt, BSplines_DTI.txt, BSplines_T1.txt, GroupWise_Huizinga_DCE.txt, GroupWise_Huizinga_DTI.txt, GroupWise_Huizinga_DCE.txt, GroupWise_Huizinga_T1.txt 
|
|
|
|
|
+-- Leeds_Patient_4128001 (each patient folder i.e. Leeds_Patient_4128001 contains folders 19, 31 and 39 which contain the acquired data for T1, DTI and DCE respectively)
|         |
|         +-- 19
|         |    |
|         |    +-- DICOM (this folder is expected to contain all DICOM files constituting the T1 acquisition, the code will split them into slices)
|         |                    
|         |
|         +-- 31
|         |    |
|         |    +-- DICOM (this folder is expected to contain all DICOM files constituting the DTI acquisition, the code will split them into slices)
|         | 
|         |                    
|         |
|         +-- 39
|              |
|              +-- DICOM (this folder is expected to contain all DICOM files constituting the DCE acquisition, the code will split them into slices)
|
|
|
|
|
+-- Leeds_Patient_4128002
|
|
|
+-- ...etc...


## Notes regarding `GroupWise_VS_MDR.py`
 
1) there is a global variable called TECHNIQUE in line 1237  +-- if it is set equal to 1 (int) it performs GroupWise model-free image registration
                                                             |
                                                             +-- if it is set equal to 2 (int) it performs MDR model-driven registration



2) there is a global variable called SEQUENCES (line 1217) which is a list of python strings and it will indicate which MR acquisitions will be motion corrected,
   the strings included in SEQUENCES can be any subset/ (or the whole set) of strings found in the global variable POOL_OF_SEQUENCES (line 1216) 


3) the global variable CORRESPONDANCE (line 1227) is a dictionary that has as key values the different acquisition mechanisms (the strings found in the variable POOL_OF_SEQUENCES (line 1216))
   
   and for each key CORRESPONDANCE contains a list of three integers.
                                                                    +-- the first shows the number of folder that corresponds to the data of the acquired sequence 
                                                                    |   (i.e. the current implementation follows the convention: T1 -> 19, DTI -> 31, DCE -> 39)
                                                                    |
                                                                    +-- the second shows the number of files for each acquired sequence
                                                                    | 
                                                                    +-- the third shows the number of slices for each acquired sequence


4) the naming correspondance between the patient folders (i.e. Leeds_Patient_001) and the folders for each patient in the AIFs folder (i.e. patient_001) is based on the last three digits of the string (i.e. 001) which are expected to be unique for each patient)

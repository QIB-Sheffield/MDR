
				                 		------------------------------------                                                                  
                               				       |                                    |
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|   Schematic of folders structure   |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                               |                                    |
                                                                ------------------------------------


Patient (-> parent folder inside it are stored all the patient data that will be processed, each patient in its own folder,
            inside the parent folder are also stored the folders: Elastix_Parameters_Files and AIFs as explained by the schematic below,
            the full path for this folder is given in line 1234 as value for the global variable: DATA_PATH)
|
|
|-> AIFs (this folder is needed for DCE data and contains the aif for each patient in a .txt form)
|    |
|    |-> patient_001
|    |          |
|    |          |-> AIF__2CFiltration__Curve.txt
|    |
|    |-> patient_002
|               |
|               |-> AIF__2CFiltration__Curve.txt
|      
|     ...etc...
|
|
|
|
|-> Elastix_Parameters_Files (the full path for this folder is given in line 1297 in the argument parameters_file=...)
|    |
|    |-> BSplines_DCE.txt, BSplines_DTI.txt, BSplines_T1.txt, GroupWise_Huizinga_DCE.txt, GroupWise_Huizinga_DTI.txt, GroupWise_Huizinga_DCE.txt, GroupWise_Huizinga_T1.txt 
|
|
|
|
|
| -> Leeds_Patient_4128001 (each patient folder i.e. Leeds_Patient_4128001 contains folders 19, 31 and 39 which contain the acquired data for T1, DTI and DCE respectively)
|         |
|         |-> 19
|         |    |
|         |    |-> DICOM (this folder is expected to contain all DICOM files constituting the T1 acquisition, the code will split them into slices)
|         |                    
|         |
|         |-> 31
|         |    |
|         |    |-> DICOM (this folder is expected to contain all DICOM files constituting the DTI acquisition, the code will split them into slices)
|         | 
|         |                    
|         |
|         |-> 39
|              |
|              |-> DICOM (this folder is expected to contain all DICOM files constituting the DCE acquisition, the code will split them into slices)
|
|
|
|
|
| -> Leeds_Patient_4128002
|
|
|
|    ...etc...
                                                                               -------
                                                                              |       |                                                                        
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~| Notes |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                              |       |
                                                                               -------

 
1) there is a global variable called TECHNIQUE in line 1236  -> if it is set equal to 1 (int) it performs GroupWise model-free image registration
                                                               |
                                                                -> if it is set equal to 2 (int) it performs MDR model-driven registration



2) there is a global variable called SEQUENCES (line 1217) which is a list of python strings and it will indicate which MR acquisitions will be motion corrected,
   the strings included in SEQUENCES can be any subset/ (or the whole set) of strings found in the global variable POOL_OF_SEQUENCES (line 1216) 


3) the global variable CORRESPONDANCE (line 1227) is a dictionary that has as key values the different acquisition mechanisms (the strings found in the variable                                                                                                                                             POOL_OF_SEQUENCES (line 1216))
   
   and for each key CORRESPONDANCE contains a list of three integers-> the first shows the number of folder that corresponds to the data of the acquired sequence 
                                                                    |   (i.e. the current implementation follows the convention: T1 -> 19, DTI -> 31, DCE -> 39)
                                                                    |
                                                                     -> the second shows the number of files for each acquired sequence
                                                                    | 
                                                                     -> the third shows the number of slices for each acquired sequence


4) the naming correspondance between the patient folders (i.e. Leeds_Patient_001) and the folders for each patient in the AIFs folder (i.e. patient_001) is based on the last three digits of the string (i.e. 001) which are expected to be unique for each patient)

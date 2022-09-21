# Compute N2 from GO-SHIP EasyOcean reported data
I downloaded the GO-SHIP EasyOcean reported data as of September 2021 (https://github.com/kkats/GO-SHIP-Easy-Ocean) <br />
1. computed N2 from the *reported data* <br />
2. interpolated onto standard grids following GO-SHIP EasyOcean *gridded data* procedures <br />
<br />
To do that, I <br />
1. clone the GO-SHIP EasyOcean repository to my local work folder: git clone https://github.com/kkats/GO-SHIP-Easy-Ocean ../my_work_folder <br />
2. replace grid_data_pressure.m in ../my_work_folder/ with the modified file in my repository: cp grid_data_pressure.m ../my_work_folder/ <br />
3. copy STrun.m to ../my_work_folder/: cp STrun.m ../my_work_folder/ <br />
4. Line 13 and 14 in STrun.m change the directory to the directory where you store the *reported data* and *gridded data*, respectively. <br />
5. run STrun.m with MATLAB and you will get *gridded data* that contains four estimates of N2, the paper use CTDN2_f4, see paper for details.

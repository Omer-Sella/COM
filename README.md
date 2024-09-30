# Channel Operating Margin (COM)
# Table of contents
1. [What this repository is about](#about)
2. [Structure](#structure)
3. [License](#license)


## What this repository is about<a name="about"></a>
COM, or, Channel Operating Margin, is part of the 802.3 Ethernet physical layer specification.\
It is an algorithm, or a calculation, which is used to determine whether an electrical channel (cable, trace, package) is compliant, and what is its expected performance.\
This repository is my first step in providing the (802.3) community automation tools, and later on porting COM to python.
## Structure<a name="structure"></a>
The folder /sourceFiles contains MATLAB files downloaded from 802.3 pages. In these files, all functions appear in one .m file\
The individual functions were separated (by me in this case) and can be found under /individualFunctions.\
Channel contributions (also downloaded from 802.3 pages) are found as **submodules** under /channels\
The majority of the work is currently under /tests, where I created 
#### 1. A "test setup" tree that reflects the directory tree in /channels, where the leaf directories contain .json files that define setups (thru, NEXT, FEXT and configuration sheet) to be used with COM.
#### 2. A results tree, again, similar to the channels tree, that would contain the results of running COM using the .json files and consequently the channels.
## License<a name="license"></a>
The COM MATLAB and Excel files in this projects were publicly available in 802.3dj tools under https://www.ieee802.org/3/dj/public/tools/index.html and other 802.3 project websites.\
No ownership or credit is claimed here, and this repository is just for contribution. Any usage, distribution, re-distribution is contingent on whatever updated license is provided in the respective tools pages - so please check there.\
The python files in this project were created by me, as a contribution to development of the standard.

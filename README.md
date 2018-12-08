# CMSSW
CMS - Software C++ and Python Scripts
Install the appropriate cms release
For 53X:
cmsrel CMSSW_5_3_32
cd CMSSW_5_3_32/src
cmsenv
git clone git@github.com:stjohnso98/CMSSW.git           # not including any corrections right now
scram b -rj16

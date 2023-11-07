#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# diff_exocam.py
#
# Author Eric Wolf
# November 2023
#
# Purpose:  diffs beteween a $CASE and the ExoCAM source code.
#           with --case2 $CASE2 option, diffs between two CASEs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



import os
import sys
import subprocess
import argparse

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("~~~~~~~~~~~~~~~~~~~ difference ExoCAM ~~~~~~~~~~~~~~~~~~~")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

                                                                    
parser = argparse.ArgumentParser()
# user required inputs
parser.add_argument('case'      , type=str,   nargs=1, default=' ', \
                     help='set case name')
parser.add_argument('--case2'      , type=str,   help='set case2 name')
args = parser.parse_args()

case      = str(args.case[0])



# read in system directories specified in sys.in 
with open('sys.in', 'r') as f:
        machine     = f.readline()
        caseroot    = f.readline()
        ec_loc      = f.readline()
        ert_loc     = f.readline()
        cesm_loc    = f.readline()

machine = machine.rstrip("\n")

########################################################################
# Set directories 
#

# Scripts directory
scriptsdir = cesm_loc + '/cesm1_2_1/scripts'
scriptsdir = ''.join(scriptsdir.split())

# CASE directories
case_dir           = caseroot + '/' + case
case_dir           = ''.join(case_dir.split())
case_cam_dir        = caseroot + "/" + case + "/SourceMods/src.cam"
case_cam_dir        = ''.join(case_cam_dir.split())
case_clm_dir        = caseroot + "/" + case + "/SourceMods/src.clm"
case_clm_dir        = ''.join(case_clm_dir.split())
case_cice_dir       = caseroot + "/" + case + "/SourceMods/src.cice"
case_cice_dir       = ''.join(case_cice_dir.split())
case_drv_dir        = caseroot + "/" + case + "/SourceMods/src.drv"
case_drv_dir        = ''.join(case_drv_dir.split())
case_share_dir      = caseroot + "/" + case + "/SourceMods/src.share"
case_share_dir      = ''.join(case_share_dir.split())

# CASE2 directories
if args.case2 is not None:
    case2_dir           = caseroot + '/' + args.case2
    case2_dir           = ''.join(case2_dir.split())
    case2_cam_dir        = caseroot + "/" + args.case2 + "/SourceMods/src.cam"
    case2_cam_dir        = ''.join(case2_cam_dir.split())
    case2_clm_dir        = caseroot + "/" + args.case2 + "/SourceMods/src.clm"
    case2_clm_dir        = ''.join(case2_clm_dir.split())
    case2_cice_dir       = caseroot + "/" + args.case2 + "/SourceMods/src.cice"
    case2_cice_dir       = ''.join(case2_cice_dir.split())
    case2_drv_dir        = caseroot + "/" + args.case2 + "/SourceMods/src.drv"
    case2_drv_dir        = ''.join(case2_drv_dir.split())
    case2_share_dir      = caseroot + "/" + args.case2 + "/SourceMods/src.share"
    case2_share_dir      = ''.join(case2_share_dir.split())
    
# ExoCAM directories
exo_cam_dir     = ec_loc + "/ExoCAM/cesm1.2.1/configs/cam_aqua_fv/SourceMods/src.cam"
exo_cam_dir     = ''.join(exo_cam_dir.split())
exo_clm_dir     = ec_loc + "/ExoCAM/cesm1.2.1/configs/cam_land_fv/SourceMods/src.clm"
exo_clm_dir     = ''.join(exo_clm_dir.split())
exo_cice_dir    = ec_loc + "/ExoCAM/cesm1.2.1/configs/cam_aqua_fv/SourceMods/src.cice"
exo_cice_dir    = ''.join(exo_cice_dir.split())
exo_drv_dir     = ec_loc + "/ExoCAM/cesm1.2.1/configs/cam_aqua_fv/SourceMods/src.drv"
exo_drv_dir     = ''.join(exo_drv_dir.split())
exo_share_dir   = ec_loc + "/ExoCAM/cesm1.2.1/configs/cam_aqua_fv/SourceMods/src.share"
exo_share_dir   = ''.join(exo_share_dir.split())

# ExoCAM subdirectories
#fv_dir         = exo_cam_dir + "/src.cam.fv"
#fv_dir         = ''.join(fv_dir.split())
#carma_dir      = exo_cam_dir + "src.carma"
#carma_dir      = ''.join(carma_dir.split())

# ExoRT radiation directories
# still not quite sure what to do with radiation files
# for now copies are in ExoCAM
# but my goal here is to have zero redunancy of files, except where absolutely neccessary
#n28archean_dir     = ert_loc + "/ExoRT/3dmodels/src.cam.n28archean"
#n28archean_dir     = ''.join(n28archean_dir.split())
#n42h2o_dir         = ert_loc + "/ExoRT/3dmodels/src.cam.n42h2o"
#n42h2o_dir         = ''.join(n42h2o_dir.split())
#n68equiv_dir       = ert_loc + "/ExoRT/3dmodels/src.cam.n68equiv"
#n68equiv_dir       = ''.join(n68equiv_dir.split())
#n68equiv_haze_dir  = ert_loc + "/ExoRT/3dmodels/src.cam.n68equiv_haze"
#n68equiv_haze_dir  = ''.join(n68equiv_haze_dir.split())

# namelist directory
#namelist_dir = ec_loc + "/ExoCAM2/namelists"

########################################################################
#

print("differencing files...")
if args.case2 is not None:
    print(case, " vs. ", args.case2)
    comp_cam_dir   = case2_cam_dir
    comp_share_dir = case2_share_dir
    comp_drv_dir   = case2_drv_dir
    comp_clm_dir   = case2_clm_dir
    comp_cice_dir  = case2_cice_dir
else:
    print(case, " vs. ExoCAM Source")
    comp_cam_dir   = exo_cam_dir
    comp_share_dir = exo_share_dir
    comp_drv_dir   = exo_drv_dir
    comp_clm_dir   = exo_clm_dir
    comp_cice_dir  = exo_cice_dir

print(" ")
print("<<< <<< <<< <<< <<< <<< <<< src.cam >>> >>> >>> >>> >>> >>> >>>")
if os.path.exists(case_cam_dir) and os.path.isdir(case_cam_dir):
    for filename in os.listdir(case_cam_dir):
        print(filename)                
        print("======================================================")
        f = ['diff', case_cam_dir + '/' + filename, comp_cam_dir + '/' + filename]
        subprocess.run(f)

print(" ")
print("<<< <<< <<< <<< <<< <<< <<< src.share >>> >>> >>> >>> >>> >>> >>>")
if os.path.exists(case_share_dir) and os.path.isdir(case_share_dir):
    for filename in os.listdir(case_share_dir):
        print(filename)                
        print("======================================================")
        f = ['diff', case_share_dir + '/' + filename, comp_share_dir + '/' + filename]
        subprocess.run(f)

print(" ")
print("<<< <<< <<< <<< <<< <<< <<< src.drv >>> >>> >>> >>> >>> >>> >>>")
if os.path.exists(case_drv_dir) and os.path.isdir(case_drv_dir):
    for filename in os.listdir(case_drv_dir):
        print(filename)                
        print("======================================================")
        f = ['diff', case_drv_dir + '/' + filename, comp_drv_dir + '/' + filename]
        subprocess.run(f)

print(" ")
print("<<< <<< <<< <<< <<< <<< <<< src.clm >>> >>> >>> >>> >>> >>> >>>")
if os.path.exists(case_clm_dir) and os.path.isdir(case_clm_dir):
    for filename in os.listdir(case_clm_dir):
        print(filename)                
        print("======================================================")
        f = ['diff', case_clm_dir + '/' + filename, comp_clm_dir + '/' + filename]
        subprocess.run(f)
else:
    print("No land model in $CASE")

print(" ")
print("<<< <<< <<< <<< <<< <<< <<< src.cice >>> >>> >>> >>> >>> >>> >>>")
if os.path.exists(case_cice_dir) and os.path.isdir(case_cice_dir):
    for filename in os.listdir(case_cice_dir):
        print(filename)                
        print("======================================================")
        f = ['diff', case_cice_dir + '/' + filename, comp_cice_dir + '/' + filename]
        subprocess.run(f)
else:
    print("No ice model in $CASE")


print("... finished differencing")
sys.exit



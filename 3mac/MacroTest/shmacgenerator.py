totalfiles = 10
foldername_mac = "MacroTest"
foldername_root = "MacroTest_output"
dir_5rootMacros = "/Users/trevorg04/G4ChargeExchange/5rootMacros"
dir_CE = "/Users/trevorg04/G4ChargeExchange/build"

filenameSH = f"run.sh"
with open(filenameSH, "a") as f:
    f.write(f"""
source ~/.zshrc
cd {dir_CE}      
""")
    
def var(x):
    return (x+1)*50+1

for i in range(totalfiles):
    variable = var(i)
    filename = f"run_{i}.mac"
    with open(filename, "w") as f:
        f.write(f"""
/run/initialize

#Detector
/detector/nCols {variable}
/detector/nRows {variable}
/detector/nCols_detector {variable}
/detector/nRows_detector {variable}
#/detector/nBlocks
#/detector/randPosDetector
#/detector/randPosScaling
/run/reinitializeGeometry

#Paritcle Position
#/randGun/randPAngle
#/randGun/randPos

#Folder
/output/folder {foldername_root}/{variable}

#Polarization
/Polarization/PPIndex 7
/Polarization/PrintPPIndex

/run/beamOn 100000
""")

for i in range(totalfiles):
    variable = var(i)
    with open(filenameSH, "a") as f:
        f.write(f"""./sim {foldername_mac}/run_{i}.mac   
""")
        
for i in range(totalfiles):
    variable = var(i)
    with open(filenameSH, "a") as f:
        f.write(f"""cp -r {dir_5rootMacros}
    {dir_CE}/root/{foldername_root}/{variable}
""")

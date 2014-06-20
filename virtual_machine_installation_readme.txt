A detailed documentation on how to USE CoverageCheck is at https://github.com/aweller/CoverageCheck
This document only explains how to SETUP CoverageCeck as a Virtual Machine on a Windows 7 PC. 

1. Download the CoverageCheck Virtual Machine Image from http://www.oxford-translational-molecular-diagnostics.org.uk/

2. Download and install VirtualBox from https://www.virtualbox.org/wiki/Downloads

3. Start VirtualBox. Go to File -> Import Appliance. Select the downloaded Virtual Machine Image and click import. 

4. CoverageCheck should now show up as a Virtual Machine (VM) on the left panel of VirtualBox. 
Select the VM, then go to Settings -> Shared Folders  and click the "blue folder with a plus" button on the right side of setup a new shared folder.
Select a folder where you want to keep your bams. This folder will later be accesible to CoverageCheck. 
Set the folder name to "share". Tick the "auto-mount" checkbox.

5. You're now ready to start the VM by selectin "Start".

6. Once the VM has booted up, you should see a desktop with several launcher icons.

- "CoverageCheck settings" lets you modify e.g. the default coverage cutoff or the plotting options. Don't forget to save after changing something!

- "Update CoverageCheck" will download the latest version from GitHub. This will also recreate your settings file. 

- "Activate shared folder" will allow you to see the shared folder from within the Linux VM. This can of course only work if you correctly setup the shared folder 	BEFORE starting the VM (see item 4.)

- "START COVERAGECHECK" ... guess what this might do? :) 

7. When you're done with analyzing the coverage, move your mouse to the bottom of the screen so the VirtualBox menu pops up.
Select "Machine" -> "Close" -> "Save the machine state". This way the VM just "freezes" instead of being shut down, so you don't need to wait or click the first icon once you return.


 

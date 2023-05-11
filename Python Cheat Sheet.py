# 11 May 2023
# Ryan Schlimme

# This scratch document provides quick tricks and code to aid in running python files through the terminal


##### TERMINAL Commands #####
ls 			# lists elements of a directory/folder currently running
cd "name"		# changes directories
cd .. 		# goes up one level in the path 
pwd 			# displays current path
clear 		# erases display of console (preserves memory of past commands)
python 		# enters python
exit()		# exits python
python "path" 	# executes .py file in the listed path (copy path w/quotes by Ctrl+Shift+C)


##### Pip Package Management #####
pip list 				# lists all installed packages
pip install --upgrade pip 	# upgrades pip
pip install name 			# installs package called "name"
pip uninstall name 		# uninstalls package called "name"


##### TDMS/Brownian Functions #####
f_name = r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\iter_0.tdms" 	# create a variable pointing to file (change Ryan Schlimme to ryans)
import sys													# allows access to path.append
sys.path.append(r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS 								# pull function from time_series module
C = CollectionTDMS(f_name) 										# class to apply various functions to time series data
C.set_collection("Y") 										 	# select channel Y

##### Classes #####
print(dir(class name)) 	# prints all methods of a class

##### Plotting #####
import matplotlib.pyplot as plt 	# imports necessary modules
plt.loglog(x,y) 				# plots on a logx,logy scale
plt.show					# display the plot


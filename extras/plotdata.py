import os
from numpy.random import rand

try:
	import matplotlib.pyplot as plt
except:
	var=raw_input("It looks like you don't have matplotlib installed, would you like to install it [Y/N]?")
	inc=0
	while True:
		if (var.upper()=="Y"):
			print "Great, attempting to install python-matplotlib."
			os.system('sudo apt-get install -y python-matplotlib')
			break
		elif (var.upper()=="N"):
			print "Ok."
			quit()
		else:
			var=raw_input("Y or N please:")
			inc+=1
		if (inc==10):
			print "Error entering data, exiting."
			quit()

print """please cite:
  @Article{Hunter:2007,
  Author    = {Hunter, J. D.},
  Title     = {Matplotlib: A 2D graphics environment},
  Journal   = {Computing In Science \& Engineering},
  Volume    = {9},
  Number    = {3},
  Pages     = {90--95},
  abstract  = {Matplotlib is a 2D graphics package used for Python
  for application development, interactive scripting, and
  publication-quality image generation across user
  interfaces and operating systems.},
  publisher = {IEEE COMPUTER SOC},
  year      = 2007"""

for color in ['red', 'green', 'blue']:
    n = 750
    x, y = rand(2, n)
    scale = 200.0 * rand(n)
    plt.scatter(x, y, c=color, s=scale, label=color,
                alpha=0.3, edgecolors='none')

plt.legend()
plt.grid(True)

plt.show()

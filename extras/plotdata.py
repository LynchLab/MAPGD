try:
	import matplotlib
except:
	try:
		import os
	except:
		print "It looks like you don't have matplotlib installed, you need a copy to run this script."
		
	var=raw_input("It looks like you don't have matplotlib installed, would you like to install it [Y/N]?")
	inc=0
	while True:
		if (var.upper()=="Y"):
			print "YAH!"
			os.system('sudo apt-get install -y python-matplotlib')
			break
		elif (var.upper()=="N"):
			print "*sigh* fine."
			quit()
		else:
			var=raw_input("Y or N please:")
			inc+=1
		if (inc==10):
			print "What are you doing, typing with your feet?"
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


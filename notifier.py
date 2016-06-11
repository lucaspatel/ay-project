#!
'''
@author: lucas
'''

from subprocess import Popen, PIPE

def promptSuccess(title, message):
	prompt(title, message, "Glass")

def promptError(title, message):
	prompt(title, message, "Basso")

def prompt(title, message, sound):
	script = "osascript -e 'display notification \""+ message + "\" with title \"" + title + "\" sound name \"" + sound + "\"'"
	p = Popen(script, shell=True,stdin=PIPE, stdout=PIPE, stderr=PIPE)
	stdout, stderr = p.communicate(script)

def main():
	prompt("Test title!","Test message!","Ping")

if __name__ == "__main__":
	main()
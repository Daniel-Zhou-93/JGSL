JAVA_SRC = src
PACKAGE = class/JGSL
SUBPACKAGES_CLASS_PATH = class
JC = javac
JFLAGS = -d

all: package 

package:
	$(MAKE) -C src
	$(MAKE) -C test
	
clean:
	$(MAKE) clean -C test
	$(MAKE) clean -C src


JAVA_HOME = ..

EXAMPLE_CLASS = JGSLExamples.class
EXAMPLE_RUN = JGSL_examples_result.txt
LU_CLASS = JGSLExamplesLU.class
LU_RUN = JGSL_examples_LU_result.txt
CHOL_CLASS = JGSLExamplesChol.class
CHOL_RUN = JGSL_examples_Chol_result.txt
EIG_CLASS = JGSLExamplesEig.class
EIG_RUN = JGSL_examples_Eig_result.txt
QR_CLASS = JGSLExamplesQR.class
QR_RUN = JGSL_examples_QR_result.txt
SVD_CLASS = JGSLExamplesSVD.class
SVD_RUN = JGSL_examples_SVD_result.txt

all: $(EXAMPLE_CLASS) $(EXAMPLE_RUN) $(LU_CLASS) $(LU_RUN) $(CHOL_CLASS) $(CHOL_RUN) $(EIG_CLASS) $(EIG_RUN) $(QR_CLASS) $(QR_RUN) $(SVD_CLASS) $(SVD_RUN)

$(EXAMPLE_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JGSLExamples.java 

$(EXAMPLE_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JGSLExamples > JGSL_examples_result.txt

$(LU_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JGSLExamplesLU.java

$(LU_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JGSLExamplesLU > JGSL_examples_LU_result.txt

$(CHOL_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JGSLExamplesChol.java

$(CHOL_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JGSLExamplesChol > JGSL_examples_Chol_result.txt

$(EIG_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JGSLExamplesEig.java

$(EIG_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JGSLExamplesEig > JGSL_examples_Eig_result.txt

$(QR_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JGSLExamplesQR.java

$(QR_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JGSLExamplesQR > JGSL_examples_QR_result.txt

$(SVD_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JGSLExamplesSVD.java

$(SVD_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JGSLExamplesSVD > JGSL_examples_SVD_result.txt

clean: 
	rm -f *.txt *.class *.out *.serial

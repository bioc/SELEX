all: xjc compile jar install

xjc:
	xjc selex.xsd -p config

compile:
	javac `find . -name '*.java'`

jar:
	jar cvf selex.jar `find . -name '*.class'`

install:
	cp -pv selex.jar ../../inst/java/selex.jar

clean:
	rm -rf config
	rm -f `find . -name '*.class'`

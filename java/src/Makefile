# run with Java 8 to target 8 and up
all: xjc compile jar install clean

xjc:
	xjc selex.xsd -p config

compile:
	mkdir -p jaxb; cd jaxb; \
	wget https://repo1.maven.org/maven2/com/sun/xml/bind/jaxb-ri/2.3.1/jaxb-ri-2.3.1.zip; \
	unzip jaxb-ri-2.3.1.zip; \
	mv jaxb-ri/mod/*.jar jaxb-ri/bin/xjc.sh jaxb-ri/*License* .; \
	rm -rf jaxb-ri jaxb-ri-2.3.1.zip
	export CLASSPATH="jaxb/jaxb-api.jar:jaxb/jaxb-runtime.jar:jaxb/javax.activation-api.jar"; \
	javac `find . -name '*.java'`

jar:
	jar cvf selex.jar `find . -name '*.class'`
	mkdir tmp
	cd tmp; unzip -uo ../selex.jar; \
	unzip -uo ../jaxb/jaxb-api.jar; \
	unzip -uo ../jaxb/javax.activation-api.jar; \
	unzip -uo ../jaxb/jaxb-runtime.jar; \
	unzip -uo ../jaxb/istack-commons-runtime.jar
	jar -cvf selex.jar -C tmp .
	rm -rf tmp

install:
	cp -pv selex.jar ../../inst/java/selex.jar
	rm -rf jaxb/*.jar

clean:
	rm -f `find . -name '*.class'`

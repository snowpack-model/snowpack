#!/bin/sh
#This is just the necessary command line to run a jar file

#launcher for alpine3d data viewer
#java -Xms128m -Xmx800m -jar view.jar

ROOT=`dirname $0`
JAR="view.jar"
EXTRA="-Xms128m -Xmx800m"

if [ -f ${ROOT}/${JAR} ]
then
	java ${EXTRA} -jar ${ROOT}/${JAR}
else if [ -d ${ROOT}/bin ]
	then
		java ${EXTRA} -jar ${ROOT}/bin/${JAR}
	else
		java ${EXTRA} -jar ${ROOT}/tools/Interface/${JAR}
	fi
fi

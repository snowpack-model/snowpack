#!/bin/sh
#This is just the necessary command line to run a jar file
JAR="inishell.jar"
OS=`uname`
ROOT=`dirname $0`

if [ "${OS}" = "Darwin" ]
then
	EXTRA="-Xdock:name=INIshell"
fi

if [ -f ${ROOT}/${JAR} ]
then
	java ${EXTRA} -jar ${ROOT}/${JAR}
else if [ -d ${ROOT}/bin ]
	then
		java ${EXTRA} -jar ${ROOT}/bin/${JAR}
	else
		java ${EXTRA} -jar ${ROOT}/dist/${JAR}
	fi
fi

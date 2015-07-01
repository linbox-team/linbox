#!/bin/bash - 

# Copyright(c) 2011 LinBox
# Written by BB <bboyer@imag.fr>
# see COPYING for more details.


STABLE_LINBOX=1.3.0
STABLE_VAR="true"
DONE="\033[0;36m done !\033[0m"

#########
#  die  #
#########

die() {
	echo -ne "\n\033[1;31m * \033[0mfailed"
	if [[ -n $1 ]] ; then
		echo " ($1)"
	else 
		echo "."
	fi
	exit -1 ;  
}

cool() {
	echo -e $DONE| tee -a auto-install.log
}

############
#  parser  # 
############


for i in "$@" ; do
	case "$i" in
	"--stable")
		if [ "x$STABLE_VAR" = "xfalse" ] ; then echo "stable or not ?";              exit -1 ;  fi
		STABLE_VAR=true;
		;;
	*)
		if [[ ! "$i" =~ --.*=.+ ]] ; then
			echo "bad switch : $i"
			help ;
			exit -1 ;
		fi
		# options (now we can cut)
		QUI="`echo $i | cut -d'=' -f1`"
		QUOI="`echo $i | cut -d'=' -f2`"
		# echo "$QUI = $QUOI"
		case "$QUI" in
			"--stable")
				OK=2
				[[ "$QUOI" = "yes" ]]  &&  OK=1  
				[[ "$QUOI" = "no"  ]]  &&  OK=0   
				if [[ "$OK" = "2" ]] ; then 
					echo "stable=[yes/no] !" ; exit -1 ; 
				fi
				[[ "OK" = "1" ]] && STABLE_VAR="true" || STABLE_VAR="false"

				;;
		esac
		;;
	esac
done


### Fetch ###

echo -en "${BEG}fecthing LinBox..."| tee -a auto-install.log
if [ "$STABLE_VAR" = "true" ]; then
	if [ -f linbox-$STABLE_LINBOX.tar.gz ] ; then
		echo "already there"
		echo -ne "${BEG}fetching md5sum" ; 
		[ -f linbox-$STABLE_LINBOX.tar.gz.md5sum ] && rm linbox-$STABLE_LINBOX.tar.gz.md5sum ;
		wget http://www.linalg.org/linbox-$STABLE_LINBOX.tar.gz.md5sum >/dev/null 2>&1 || die
		[ -f linbox-$STABLE_LINBOX.tar.gz.md5sum  ] || die
		cool
		echo -ne "${BEG}"
		md5sum -c linbox-$STABLE_LINBOX.tar.gz.md5sum || die
	else
		wget http://www.linalg.org/linbox-$STABLE_LINBOX.tar.gz >/dev/null 2>&1 || die
		[ -f linbox-$STABLE_LINBOX.tar.gz ] || die
		echo -ne "${BEG}fetching md5sum" ; 
		wget http://www.linalg.org/linbox-$STABLE_LINBOX.tar.gz.md5sum >/dev/null 2>&1 || die
		cool
		echo -ne "${BEG}"
		md5sum -c linbox-$STABLE_LINBOX.tar.gz.md5sum || die
	fi
	cool
else
	OK=0 ;
	svn co svn://linalg.org/linalg/trunk/linbox 2>&1 >/dev/null && OK=1 
	[ "$OK" = "1" ] &&  cool  || die 
fi

### Extract ###

OK=0
if [ "$STABLE_VAR" = "true" ]; then
	echo -en "${BEG}extracting LinBox..."| tee -a auto-install.log
	tar xzf linbox-$STABLE_LINBOX.tar.gz  && OK=1
	[ "$OK" = "1" ] &&  cool   || die 
fi

### Fetch the other packages and install the whole lot ...
if [ "$STABLE_VAR" = "true" ]; then
	cd linbox-$STABLE_LINBOX || die
else
	cd linbox/ || die
fi

echo -e "`pwd`/auto-install.sh $*" >> auto-install.log
./auto-install.sh $*

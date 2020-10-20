#!/bin/bash
#executing this tars up the code directory.

OS=$(uname)
case "$OS" in
	"Darwin")
	    echo "OS: Darwin, MacOS"
	    scriptdir="$(pwd "$(dirname "${BASH_SOURCE[0]}")")/"
	    codedir=$(sed 's|tools/||' <<<"$scriptdir")
	    tarname=$(sed 's|/$||;s|.*/D||' <<<"$codedir")".tar.gz"
	    ;;
	"Linux")
	    echo "OS: Linux"
	    scriptdir="$(readlink -f "$(dirname "${BASH_SOURCE[0]}")")/"
	    codedir=$(sed 's|tools/||' <<<"$scriptdir")
	    tarname=$(sed 's|/$||;s|.*/D\?||' <<<"$codedir")".tar.gz"
	    ;;
	*)
	    echo "Operating system ($OS) could not be defined"	    
	    exit 1
esac

#scriptdir="$(readlink -f "$(dirname "${BASH_SOURCE[0]}")")/"
#codedir=$(sed 's|tools/||' <<<"$scriptdir")
codeend=$(sed 's|/.*/\(.*\)/|\1|' <<<"$codedir")
tardir=$(sed 's|/$||;s|\(.*/\).*|\1|' <<<"$codedir")
tarpath="$tardir$tarname"

#cd $codedir

#echo "tar -czvf $tarpath -C $tardir $codeend"
#tar --exclude-backups -czvf $tarpath -C $tardir $codeend #$scriptdir
tar --exclude="*#" --exclude="*.#*" --exclude="*~" --exclude="TODO" --exclude="OLD_*" --exclude="*.aux" --exclude="*.log" --exclude="*.out" --exclude="*.toc" --exclude="*.txt" --exclude="${codeend}/bin" --exclude="*.tar.gz" --exclude="${codeend}/src/archive" --exclude="${codeend}/obj" --exclude="${codeend}/.git" --exclude="*.gitignore" -czvf $tarpath -C $tardir $codeend
mv $tarpath ../.
echo "made tar file $tarname"

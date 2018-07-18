#!/usr/bin/sh

if test -r "lb-maple.mpl.bak"; then
    mv lb-maple.mpl.bak lb-maple.mpl
fi

sed -e "s|lbpathvalue|lbpath:=\"$1\";|" -i.bak lb-maple.mpl  
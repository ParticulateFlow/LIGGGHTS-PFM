# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# arg1 = file, arg2 = file it depends on

action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    if (test -z "$2" || test -e ../$2) then
      cp $1 ..
      if (test $mode = 2) then
        echo "  updating src/$1"
      fi
    fi
  elif (test -n "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
    fi
  fi
}

# step 1: install partitioner

action partitioner_zoltan.h
action partitioner_zoltan.cpp

# step 2: handle cases and tasks not handled in step 1.

if (test $mode = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*ZOLTAN[^ \t]* //' ../Makefile.package
    sed -i -e 's/[^ \t]*zoltan[^ \t]* //' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-DLMP_USER_ZOLTAN |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-lzoltan |' ../Makefile.package
  fi

  # force rebuild

  touch ../atom.cpp

elif (test $mode = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*ZOLTAN[^ \t]* //' ../Makefile.package
    sed -i -e 's/[^ \t]*zoltan[^ \t]* //' ../Makefile.package
  fi

  # force rebuild of files with LMP_USER_OMP switch

  touch ../atom.cpp

fi

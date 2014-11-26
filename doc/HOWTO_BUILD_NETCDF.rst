How to build the NetCDF4 C++ development packages under Debian and Ubuntu
=========================================================================

At the moment (as of November 2014), Debian has a bug in it's NetCDF4 development pacakges in that it doesn't build the C++ bindings for NetCDF version 4.

https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=728498

The package *should* inculde the file /usr/lib/libnetcdf_c++-4.a, but it doesn't.

You can two options to get this file on your system.

Option A: Rebuild the deb package
---------------------------------

Advantages: Makes your system easier to maintain, takes care of all the dependencies for you, and tracks all of the files that get installed (so that they're easily upgradable or uninstallable). Disadvantages: Debian/Ubuntu dev tools are pretty scattered, and it's more work to get everything built correctly.

For either distro, pull down the source package:

$ apt-get source libnetcdfc++4
$ cd netcdf-4.1.3

and modify the line that defines DEB_CONFIGURE_EXTRA_FLAGS and make sure that in includes --enable-cxx and --enable-cxx-4.

$ emacs debian/rules

It should look like this:

DEB_CONFIGURE_EXTRA_FLAGS := --enable-shared --enable-pic --enable-docs-install --enable-netcdf-4  --with-hdf5=/usr --with-libcf --enable-dap --disable-dap-remote-tests --enable-cxx --enable-cxx-4

Add a line to the debian/changelog file with an incremented revision.

Then rebuild the package. How you do this depends on whether you're on Debian or Ubuntu. Debian is a little easier.

$ dpkg-buildpackage -uc -b -j16

This will take a few minutes, but if all goes well, you should have a fresh libnetcdf-dev_4.1.3-7_amd64.deb package (as well as others), one directory up (../). Install them and you should be good to go.

$ dpkg -i ../libnetcdf-*.deb

Under Ubuntu you have to work a little harder to get the dev environment set up. In Ubuntu, it's advised that you actually install a virtual environment in which to install the package. 
Follow this: http://packaging.ubuntu.com/html/getting-set-up.html and the pbuilder HOWTO: https://wiki.ubuntu.com/PbuilderHowto

Create the virtual environment, and make sure to include universe, because the required dependencies libdap-dev and libhdf5-dev are both from universe.

$ sudo pbuilder create --distribution trusty --components "main universe"
hiebert@aether:~/tmp/netcdf-4.1.3$ debuild -S
hiebert@aether:~/tmp/netcdf-4.1.3$ sudo pbuilder build ../netcdf_4.1.3-8ubuntu2.dsc

Then everything you need should be in /var/cache/pbuilder/result

Option B: Build from source tarball
-----------------------------------

Advantages: less steps; you can do it as a regular user. Disadvantage: makes your system hard to maintain because you have libraries that are not tracked and are not upgraded automatically.

To do this, get the latest source from Unidata:

https://github.com/Unidata/netcdf-cxx4/releases/tag/v4.2.1

The newest version doesn't even inlclude the --enable-cxx-4 in the configure script. I think it just has it enabled, no matter what.

Configure it to be installed in your home directory:

$ ./configure --prefix=$HOME

And build and install

$ make
$ make install


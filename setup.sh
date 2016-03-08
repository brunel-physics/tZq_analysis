export LIBCONFIG_PATH=$(find ~ /lib /usr/lib /lib64 /usr/lib64 -name libconfig++.so -print0 -quit | xargs -0 -I % dirname %);
export LIBLHAPDF_PATH=$(find ~ /lib /usr/lib /lib64 /usr/lib64 -name libLHAPDF.so -print0 -quit | xargs -0 -I % dirname %);
export INCCONFIG_PATH=$(find ~ /usr/include -name libconfig.h++ -print0 -quit | xargs -0 -I % dirname %);
export INCLHAPDF_PATH=$(find ~ /usr/include -name LHAPDF.h -print0 -quit | xargs -0 -I % dirname %);

echo "# To avoid having to run this script again, add the following lines to"
echo "# your shell's configuration file"
echo "export LIBCONFIG_PATH=$LIBCONFIG_PATH"
echo "export LIBLHAPDF_PATH=$LIBLHAPDF_PATH"
echo "export INCCONFIG_PATH=$INCCONFIG_PATH"
echo "export INCLHAPDF_PATH=$INCLHAPDF_PATH"

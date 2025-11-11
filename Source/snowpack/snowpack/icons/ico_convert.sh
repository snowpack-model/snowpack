#This is to convert several graphics files into a Windows .ico file.
#First, prepare png with transparency at the following resolutions: 16x16, 24x24, 32x32, 48x48, 64x64, 256x256
#Then, prepare the same image sizes but at 8 bit resolution without transparency, saved as bmp
#Then, use the following command to pack the files in one .ico file with a relatively new version of ImageMagic

#convert snowpack*.bmp snowpack*.png snowpack.ico

convert snowpack16.png snowpack24.png snowpack32.png snowpack48.png -alpha on -channel RGBA -depth 8 snowpack.ico

#see also http://www.iconconstructor.com/tutorials/what_is_icon/

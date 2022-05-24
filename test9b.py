# https://stackoverflow.com/questions/37783060/how-to-display-graph-and-video-file-in-a-single-frame-window-in-python

# Version 1.3
# Uso de asi._get_video_data en Modos Video y Cube
#
# - StatusBar
# - uso de Bandwidth Overload
#
#
# Starts from opencv_wx_asi290mm.py
# 2021-02-09; mod AT 2021-03-18
# Funciona con ASI290MM
# Baje ASIStudio_V1.3.2.dmg del sitio ZWO
# lo descomprimi en /Users/edisonbustos/Downloads/ASI_linux_mac_SDK_V1.16.3/
# en lib/mac hay 3 archivos que copie (sudo) en /opt/local/lib  pues alli esta libusb-10.0.0.dylib
# Esta linea indica a pycharm que hay que bajar este paquete ...import zwoasi as asi

"""
    Inputs:
    exposure time (in microseconds), gain, brightness, ROI Size and # of exposures (used only with button "Cube")
    Controls:
    Stop/Start Video toggle button, Snapshot button, Cube button and Exit button.

    Summary of operation:
    - The software always starts taking video with default values.
    - The toggle button initiates in "Stop Video" action.
    - The video must be stopped and start in order to take a new set of input parameters
    - The "Snapshot" stops the video timer, then it sets the camera with new parameters. It takes one image,
    which is displayed and saved the image a FITS file.
    - The "Cube" stops the video timer, then it sets the camera with new parameters. It takes a "Data Cube" with the
    number of exposures already set.
    - In order to resume Video, press Stop Video and Start Video.
"""

import os
import sys
from datetime import datetime

import cv2
import wx
import zwoasi as asi
from astropy.io import fits
import numpy as np


class ASIcamera:
    def __init__(self):
        dirname = os.path.abspath('/home/andrei/ASIStudio/lib/')  # Andrei's linux laptop
        # dirname = os.path.abspath('/usr/local/lib/')              # NUC pc
        # dirname = os.path.abspath('/opt/local/lib/')  # Edison's macbook

        library_file = os.path.join(dirname, 'libASICamera2.so')  # under linux
        #library_file = os.path.join(dirname, 'libASICamera2.dylib')  # under Mac OS

        if not os.path.exists(library_file):
            print('library file not found')
            sys.exit(0)
        asi.init(library_file=library_file)

        num_cameras = asi.get_num_cameras()
        if num_cameras == 0:
            print('No cameras found')
            sys.exit(0)
        txt = '%d camera(s) found.' % num_cameras
        # print(txt)

        cameras_found = asi.list_cameras()
        camera_id = 0
        self.status = txt + '  ' + 'Camera #%d: %s' % (camera_id, cameras_found[camera_id])
        # print(txt)

        self.cap = asi.Camera(camera_id)
      
        return

    def stop(self):
        try:
            self.cap.stop_video_capture()
            self.cap.stop_exposure()
            self.video = 0
        except (KeyboardInterrupt, SystemExit):
            raise
        except (asi.ZWO_Error, asi.ZWO_IOError, asi.ZWO_CaptureError) as e:
            print('camera error:{}'.format(e))

    def set(self, exposure, gain, brightness, bandwidth, x0, y0, width, height):
        # Stop all exposures
        self.stop()

        self.cap.disable_dark_subtract()
        self.cap.set_control_value(asi.ASI_BANDWIDTHOVERLOAD, bandwidth)  # 40 for video, 80 for frame
        iexp = int(1000*exposure)
        self.cap.set_control_value(asi.ASI_EXPOSURE, iexp) # in microseconds
        self.cap.default_timeout = (self.cap.get_control_value(asi.ASI_EXPOSURE)[0] / 1000) * 2 + 500
        # self.cap.set_control_value(asi.ASI_FLIP, self.flip)
        self.cap.set_control_value(asi.ASI_GAIN, gain)
        self.cap.set_control_value(asi.ASI_HIGH_SPEED_MODE, 1) # use 1 for high-speed? 
        self.cap.set_control_value(asi.ASI_BRIGHTNESS, brightness)
        self.cap.set_roi(start_x=x0, start_y=y0, width=width, height=height, bins=1, image_type=asi.ASI_IMG_RAW16)

    def check(self, msg):
        print(msg)
        print('Exposure: {}'.format(self.cap.get_control_value(asi.ASI_EXPOSURE)[0]))
        print('Gain: {}'.format(self.cap.get_control_value(asi.ASI_GAIN)[0]))
        print('Brightness: {}'.format(self.cap.get_control_value(asi.ASI_BRIGHTNESS)[0]))
        print('Bandwidth Overload: {}'.format(self.cap.get_control_value(asi.ASI_BANDWIDTHOVERLOAD)[0]))
        print('High Speed: {}'.format(self.cap.get_control_value(asi.ASI_HIGH_SPEED_MODE)[0]))
        print('Binning:{}'.format(self.cap.get_bin()))
        print('ROI(x0,y0)(w,h):({},{})({},{})'.format(self.cap.get_roi()[0], self.cap.get_roi()[1],
                                                      self.cap.get_roi()[2], self.cap.get_roi()[3]))


class VideoPanel(wx.Panel):

    def __init__(self, parent, size):
        wx.Panel.__init__(self, parent, -1, size=size)
        self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.parent = parent
        self.SetDoubleBuffered(True)
        self.x0 = 0
        self.y0 = 0

    def OnPaint(self, event):
        dc = wx.BufferedPaintDC(self)
        dc.Clear()
        if self.parent.bmp:
            dc.DrawBitmap(self.parent.bmp, self.x0, self.y0)


def onclickExit(e):
    wx.Exit()


class MyFrame(wx.Frame):

    def __init__(self):
        wx.Frame.__init__(self, None, title="RINGSS acq software beta")

        self.bmp = None
        self.APP_TITLE = os.path.basename(__file__)
        self.Site = 'La Serena'

        self.oldtime = datetime.utcnow()

        # Frame display defaults
        self.framewidth = 512
        self.frameheight = 512

        # Video/Snapshot & Cube Defaults
        self.exposure_time = 50 #50ms in full field
        self.gain = 200
        self.brightness = 10
        self.bandwidth = 80  # 40 or 80
        self.roi_width = 1024
        self.roi_height = 1024
        self.roi_x0 = int((1920 - self.roi_width) / 2)
        self.roi_y0 = int((1080 - self.roi_height) / 2)
        self.number = 2000
        self.max_attempts = 10
        self.savedir ='data/'
        self.star = 'none'

        self.videopPanel = VideoPanel(self, (self.framewidth, self.frameheight))

        # Controls

        #font = wx.Font(pointSize=12, family=wx.FONTFAMILY_MODERN, style=wx.FONTSTYLE_ITALIC, weight=wx.FONTWEIGHT_BOLD)
        
        st_gain = wx.StaticText(self, -1, "Gain")
        self.tctrl_gain = wx.TextCtrl(self, -1,  value=str(self.gain))
        #tctrl_gain.SetFont(font)
        #tctrl_gain.Bind(wx.EVT_TEXT, self.OnGainTyped)

        st_star = wx.StaticText(self, -1, "Star")
        self.tctrl_star = wx.TextCtrl(self, -1, value=self.star)
        #tctrl_star.Bind(wx.EVT_TEXT, self.OnStarTyped)
 
        st_number = wx.StaticText(self, -1, "# of exposures")
        self.tctrl_number = wx.TextCtrl(self, -1, str(self.number))
        #tctrl_number.Bind(wx.EVT_TEXT, self.OnnumberTyped)
       
        btn_1024= wx.Button(self, id=-1, size=(100,-1), label="1024")
        btn_1024.Bind(wx.EVT_BUTTON, self.OnClick1024)

        btn_256= wx.Button(self, id=-1, size=(100,-1),label="256")
        btn_256.Bind(wx.EVT_BUTTON, self.OnClick256)

        btn_64= wx.Button(self, id=-1,size=(100,-1), label="64")
        btn_64.Bind(wx.EVT_BUTTON, self.OnClick64)

        self.format = wx.StaticText(self, -1, "1024/50ms") # current mode
                
       
        self.tbtn_video = wx.ToggleButton(self, id=-1,size=(100,-1), label="Stop Video")
        self.tbtn_video.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggleVideo)
   
        btn_frame = wx.Button(self, id=-1, label="Snapshot")
        btn_frame.Bind(wx.EVT_BUTTON, self.onclickSnapshot)

        btn_cube = wx.Button(self, id=-1, label="Cube")
        btn_cube.Bind(wx.EVT_BUTTON, self.onclickCube)

        btn_exit = wx.Button(self, id=-1, label="Exit")
        btn_exit.Bind(wx.EVT_BUTTON, onclickExit)

        self.sb = self.CreateStatusBar(4, style=wx.STB_DEFAULT_STYLE)
        self.sb.SetStatusWidths([300, 300, 300, -1])
        self.sb.SetStatusText('Video started on entry', 1)


        hbox_gain = wx.BoxSizer(wx.HORIZONTAL)
        hbox_gain.Add(st_gain, flag=wx.RIGHT, border=103)
        hbox_gain.Add(self.tctrl_gain, proportion=1)

        hbox_star = wx.BoxSizer(wx.HORIZONTAL)
        hbox_star.Add(st_star, flag=wx.RIGHT, border=103)
        hbox_star.Add(self.tctrl_star, proportion=1)

        hbox_numb = wx.BoxSizer(wx.HORIZONTAL)
        hbox_numb.Add(st_number, flag=wx.RIGHT, border=103)
        hbox_numb.Add(self.tctrl_number, proportion=1)

         
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        hbox1.Add(self.tbtn_video, 0, wx.ALIGN_CENTER)
        hbox1.Add(btn_frame, 0, wx.ALIGN_CENTER)
        hbox1.Add(btn_cube)
        
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        hbox2.Add(btn_1024, 0, wx.ALIGN_CENTER)
        hbox2.Add(btn_256, 0, wx.ALIGN_CENTER)
        hbox2.Add(btn_64, 0, wx.ALIGN_CENTER)

        hbox_format =  wx.BoxSizer(wx.HORIZONTAL)
        hbox_format.Add(self.format, flag=wx.RIGHT, border=103)
        
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.AddSpacer(10)
        vbox.Add(hbox_gain, flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, border=10)
        vbox.Add(hbox_star, flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, border=10)
        vbox.Add(hbox_numb, flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, border=10)
        vbox.AddSpacer(10)

        vbox.Add(hbox2)
        vbox.Add(hbox_format)
        vbox.AddSpacer(30)
        vbox.Add(hbox1)
        vbox.AddSpacer(10)
        vbox.Add(btn_exit, flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM, border=20)

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.videopPanel, flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM, border=40)
        sizer.Add(vbox)

        self.SetSizer(sizer)
        # self.Fit()   # This line works under Mac OS but no with
        sizer.Fit(self)
        self.Show(True)
        self.Layout()  # This is new too

        self.camera = ASIcamera()
        #self.camera.check("Debug:")
        self.sb.SetStatusText(self.camera.status)
        self.camera.set(self.exposure_time, self.gain, self.brightness, self.bandwidth, self.roi_x0, self.roi_y0,
                        self.roi_width, self.roi_height)
        self.camera.check("After Set():")

        # After camera set, these parameter must be set to use with asi._get_video_data()
        whbi = self.camera.cap.get_roi_format()
        self.sz = whbi[0] * whbi[1] * 2
        self.shape = [whbi[1], whbi[0]]  # for later
        self.timeout = (self.camera.cap.get_control_value(asi.ASI_EXPOSURE)[0] / 1000) * 2 + 500
        # Now capture can be started
        self.camera.cap.start_video_capture()

        # timer def was moved down here just after the whole setup
        self.videotimer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.OnUpdateVideo, self.videotimer)
        self.videotimer.Start(100)     #100ms refresh time in video mode 

    def OnUpdateVideo(self, event):
        timestamp = datetime.utcnow()
        # ########Debug
        aux = timestamp - self.oldtime
        self.oldtime = timestamp
        ticktime = aux.seconds + aux.microseconds / 1e6
        ##############

        # Create buffer array
        image_buffer = bytearray(self.sz)

        try:           
            asi._get_video_data(0, self.timeout, buffer_=image_buffer)           
        except (asi.ZWO_Error, asi.ZWO_IOError, asi.ZWO_CaptureError) as e:
            self.sb.SetStatusText('camera _get_video_data error {}'.format(e), 0)
            # self.videotimer.Stop()

        image_data = np.ndarray(buffer=image_buffer, dtype=np.uint16, shape=(self.shape[0], self.shape[1], 1)) // 16
        self.FrameDisplay(image_data) # display frame
        

    def OnClick1024(self, event):
        self.roi_width, self.roi_height = 1024, 1024
        self.roi_x0 = int((1920 - self.roi_width) / 2)
        self.roi_y0 = int((1080 - self.roi_height) / 2)
        self.exposure_time = 50
        self.format.SetLabel("1024/50ms")
        
    def OnClick256(self, event):
        self.roi_width, self.roi_height = 256, 256
        self.roi_x0 = int((1920 - self.roi_width) / 2)
        self.roi_y0 = int((1080 - self.roi_height) / 2)
        self.exposure_time = 10
        self.format.SetLabel("256/10ms")
        
    def OnClick64(self, event):
        self.roi_width, self.roi_height = 64, 64
        self.roi_x0 = int((1920 - self.roi_width) / 2)
        self.roi_y0 = int((1080 - self.roi_height) / 2)
        self.exposure_time = 1 # 1ms
        self.format.SetLabel("64/1ms")

    def VideoOff(self):
       self.videotimer.Stop()
       self.camera.cap.stop_video_capture()
       self.vide0 = 0
       self.sb.SetStatusText('Video stopped', 1)
       #event.GetEventObject().SetLabel("Start Video")
       self.tbtn_video.SetLabel("Start Video")
       #self.tbtn_video.SetValue(False)
       self.sb.SetStatusText(" ", 2)
       self.sb.SetStatusText(" ", 3)

    def VideoOn(self):
       # Apply GUI settings to camera
       self.GetValues()
       self.camera.set(self.exposure_time, self.gain, self.brightness, self.bandwidth, self.roi_x0, self.roi_y0,
                            self.roi_width, self.roi_height)
       self.camera.check("Set camera for Video:")

       # After camera set, these parameter must be set to use with asi._get_video_data()
       whbi = self.camera.cap.get_roi_format()
       self.sz = whbi[0] * whbi[1] * 2
       self.shape = [whbi[1], whbi[0]]  # for later
       self.timeout = (self.camera.cap.get_control_value(asi.ASI_EXPOSURE)[0] / 1000) * 2 + 500
       # Now capture can be started
       self.camera.cap.start_video_capture()

       self.videotimer.Start(100)
       self.sb.SetStatusText('Video started', 1)
       #event.GetEventObject().SetLabel("Stop Video")
       self.tbtn_video.SetLabel("Stop Video")
       #self.tbtn_video.SetValue(True)
 
    def OnToggleVideo(self, event):
        state = event.GetEventObject().GetValue()
        #print("Video: {}",self.video)
        if state:
            self.VideoOff()
        else:
            self.VideoOn()

    def GetValues(self):
        try:
            self.gain = int(self.tctrl_gain.GetValue())
            #print("gain:{}".format(self.gain))
        except ValueError:
            print("Error setting Gain, using {}".format(self.gain))
        try:
            self.star = self.tctrl_star.GetValue()
            #print("Star: "+self.star)
        except ValueError:
            text = "Wrong value for Star"
        try:
            self.number =  int(self.tctrl_number.GetValue())
            #print("Number-exp:{}".format(self.number))
        except ValueError:
            text = "Wrong value for # of exposures"

    def FrameDisplay(self, image_data): # display frame
        imax = np.max(image_data)
        imin =  np.min(image_data)          
        self.sb.SetStatusText('Imax={} Imin={}'.format(imax, imin), 2)    
        frame16 = cv2.cvtColor(image_data, cv2.COLOR_GRAY2BGR)  # frame16.dtype es uint16 y frame16.shape es 3-ch
        frame8 = cv2.convertScaleAbs(frame16, alpha=(255 / imax))  
        
        new_width = self.framewidth
        new_height = self.frameheight
        frame = cv2.resize(frame8, (new_width, new_height), interpolation=cv2.INTER_NEAREST)
        img_buf = wx.ImageFromBuffer(frame.shape[1], frame.shape[0], frame)
        self.bmp = wx.Bitmap(img_buf)  # self.bmp = wx.BitmapFromImage(img_buf)  is  deprecated
        self.Refresh()

    def onclickSnapshot(self, e):
        if self.videotimer.IsRunning():
            self.sb.SetStatusText('Timer is running! Stop Video Timer first before take a snapshot', 3)
            #self.videotimer.Stop()
            self.VideoOff()

        # Apply GUI settings to camera
        self.GetValues()
        self.camera.set(self.exposure_time, self.gain, self.brightness, self.bandwidth, self.roi_x0, self.roi_y0,
                        self.roi_width, self.roi_height)
        self.camera.check("Set camera for SnapShot:")

        timestamp = datetime.utcnow()
        try:
            image_data = self.camera.cap.capture() // 16   # image_data.dtype is uint16
        except (asi.ZWO_Error, asi.ZWO_IOError, asi.ZWO_CaptureError) as e:
            self.sb.SetStatusText('camera capture() error: {}'.format(e), 0)
            return

        self.FrameDisplay(image_data) # display frame
  
        filename = self.savedir+timestamp.strftime('%Y-%m-%d-%H%M%S') + '_snapshot' + '.fits'
        self.SaveFits(image_data, timestamp, filename)  # it writes original data
        self.VideoOn() # resume video

    def onclickCube(self, e):
        if self.roi_width > 256:
            print("No cube in full-field mode!")
            return
        if self.videotimer.IsRunning():
            self.sb.SetStatusText('Timer is running! Stop Video Timer first before take a data cube', 3)
            #self.videotimer.Stop()
            self.VideoOff()
            
        # Apply GUI settings to camera
        self.GetValues()
        self.camera.set(self.exposure_time, self.gain, self.brightness, self.bandwidth, self.roi_x0, self.roi_y0,
                        self.roi_width, self.roi_height)
        self.camera.check("Set camera for cube:")

        N = self.number
        timestamp = [None] * N
        # Create buffer array
        images_buffer = [None] * N
        # set buffer size
        whbi = self.camera.cap.get_roi_format()
        sz = whbi[0] * whbi[1] * 2
        shape = [whbi[1], whbi[0]]  # for later
        for i in range(N):
            images_buffer[i] = bytearray(sz)

        # get timeout
        timeout = (self.camera.cap.get_control_value(asi.ASI_EXPOSURE)[0] / 1000) * 2 + 500

        # Video capture
        dropped = 1
        current_attempt = 0
        for attempt in range(self.max_attempts):
            current_attempt = attempt

            self.camera.cap.start_video_capture()

            for i in range(N):
                timestamp[i] = datetime.utcnow()

                try:
                    # this is the fastest way to get the buffer data that minimizes drop frames
                    asi._get_video_data(0, timeout, buffer_=images_buffer[i])
                except (KeyboardInterrupt, SystemExit):
                    raise
                except (asi.ZWO_Error, asi.ZWO_IOError, asi.ZWO_CaptureError) as e:
                    self.sb.SetStatusText('camera _get_video_data error {}'.format(e), 0)
                    break

                if self.camera.cap.get_dropped_frames() > 0:
                    aux = (i, N, self.camera.cap.get_dropped_frames())
                    print('cube failed at iteration %d of %d with %d dropped frames' % aux)
                    break

            dropped = self.camera.cap.get_dropped_frames()
            self.camera.cap.stop_video_capture()

            if dropped == 0:
                break

        if current_attempt == (self.max_attempts - 1):
            self.sb.SetStatusText('Video failed', 3)
            return

        self.sb.SetStatusText('data cube ends, %d dropped frames of %d at attempt %d of %d' % (dropped, N,
                                                                                               current_attempt + 1,
                                                                                               self.max_attempts), 1)
 
        # Now I will process (?) and save data
        # get images from buffer with shape defined up
        img_list = [None] * N
        for i in range(N):
            img_list[i] = np.frombuffer(images_buffer[i], dtype=np.uint16).reshape(shape) // 16

        #        self.image_segmentation(img_list[N-1])

        # we skip the first value of the timestamps. No stamps if N<5
        fps = 0
        if N > 5:
            ts = np.zeros(N - 2)
            for i in range(N - 2):
                aux = timestamp[i + 2] - timestamp[i + 1]
                ts[i] = aux.seconds + aux.microseconds / 1e6

            fps = 1 / ts.mean()
            print('fps {:.1f}'.format(fps))

        # Cast the list into a numpy array
        img_array = np.array(img_list)

        filename = self.savedir+timestamp[N - 1].strftime('%Y-%m-%d-%H%M%S') + '_cube' + '.fits'

        self.SaveFits(img_array, timestamp[N - 1], filename, nframes=N, fps=fps)
        self.VideoOn() # resume video
        
    def SaveFits(self, image, tstamp, filename, nframes=0, fps=0):
        # set fits header
        hdr = fits.Header()
        hdr['EXPOSURE'] = (str(self.camera.cap.get_control_value(asi.ASI_EXPOSURE)[0]), 'Exposure time in usec')
        hdr['GAIN'] = (str(self.camera.cap.get_control_value(asi.ASI_GAIN)[0]), 'The ratio of output / input')
        hdr['BRIGHT'] = (str(self.camera.cap.get_control_value(asi.ASI_BRIGHTNESS)[0]), 'Brightness or Offset of image')
        hdr['BANDWDTH'] = (str(self.camera.cap.get_control_value(asi.ASI_BANDWIDTHOVERLOAD)[0]), 'Bandwidth Overload')
        hdr['FLIP'] = (str(self.camera.cap.get_control_value(asi.ASI_FLIP)[0]), 'Flip mode')
        hdr['HIGHSPD'] = (str(self.camera.cap.get_control_value(asi.ASI_HIGH_SPEED_MODE)[0]), 'High speed mode')
        hdr['START_X'] = (str(self.camera.cap.get_roi()[0]), 'ROI X origin in pixels')
        hdr['START_Y'] = (str(self.camera.cap.get_roi()[1]), 'ROI Y origin in pixels')
        hdr['WIDTH'] = (str(self.camera.cap.get_roi()[2]), 'ROI width in pixels')
        hdr['HEIGHT'] = (str(self.camera.cap.get_roi()[3]), 'ROI height in pixels')
        hdr['BIN'] = (str(self.camera.cap.get_bin()), 'Binning factor x,y')
        hdr['CAMERA'] = (asi.list_cameras()[0], 'Camera Model')
        hdr['CAMTEMP'] = (str(self.camera.cap.get_control_value(asi.ASI_TEMPERATURE)[0] / 10),
                          'Camera Temperature in Celsius deg')
        hdr['COMMENT'] = (str(self.camera.cap.get_camera_property()))
        hdr['SITE'] = (self.Site, 'Site Location')
        hdr['DATE-OBS'] = (str(tstamp.strftime('%Y.%m.%dT%H:%M:%S.%f')), 'UTC start date of observation')
        hdr['SWCREATE'] = (self.APP_TITLE, 'Name of software')
        hdr['NFRAMES'] = (str(nframes), 'Data cube total frames')
        hdr['FPS'] = ('{:.1f}'.format(fps), 'Frames per second')
        hdr['STAR'] = self.star

        fits.writeto(filename, image, hdr)
        txt = filename + ' saved'
        self.sb.SetStatusText(txt, 1)

if __name__ == '__main__':
    app = wx.App(0)
    myframe = MyFrame()
    app.MainLoop()

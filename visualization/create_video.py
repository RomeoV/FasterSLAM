
import numpy as np
import cv2
import matplotlib.pyplot as plt

def save_video(frames, video_path,fr=20):
    height, width, layers = np.array(frames[0]).shape
    #framerate (frames/second)
    video= cv2.VideoWriter(video_path, cv2.VideoWriter_fourcc(*'mp4v'), 2, (width,height))
    for image in frames:
        video.write(cv2.cvtColor(np.array(image), cv2.COLOR_RGB2BGR))
    video.release()
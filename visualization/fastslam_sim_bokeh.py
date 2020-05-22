

# pip3 install bokeh (!)
import time
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pprint import *
from bokeh.io import push_notebook, show, output_notebook, export_png
from bokeh.io.export import get_screenshot_as_png
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.plotting import figure 
import shutil

from selenium import webdriver

from create_video import save_video

import os
print(os.getcwd())


with open("../src/build/robot_trace_old.json") as json_file:
    robot_estimation_trace = json.load(json_file)


with open("../src/build/ground_truth_old.json") as json_file:
    ground_truth = json.load(json_file)
    
n_iters = len(ground_truth["timesteps"])



def covariance_ellipse(cov,scale=1.0):
    #Change scale here to increase size of ellipses at landmarks
    lambda_, v = np.linalg.eig(cov)
    lambda_ = np.sqrt(lambda_)
    width=lambda_[0]*2*scale
    height=lambda_[1]*2*scale
    angle=np.arccos(v[0, 0])
    return [width, height, angle]
    


def update_observations(index, trace=True):
    update_lm_dict = {
    "x": landmark_pose_trace[index][:,0] if landmark_pose_trace[index].size>0 else [],
    "y": landmark_pose_trace[index][:,1] if landmark_pose_trace[index].size>0 else [],
    "width": landmark_ellipse_trace[index][:,0] if landmark_pose_trace[index].size>0 else [],
    "height": landmark_ellipse_trace[index][:,1] if landmark_pose_trace[index].size>0 else [],
    "angle": landmark_ellipse_trace[index][:,2] if landmark_pose_trace[index].size>0 else [],
    "degree": 180*landmark_ellipse_trace[index][:,2]/np.pi if landmark_pose_trace[index].size>0 else [],
    "color" : known_landmark_colors[index],
    "idf": idf[:len(landmark_pose_trace[index])],
    }
    estimate = robot_estimation_trace["timesteps"][index]
    known_landmark_estimates_cds.data = update_lm_dict
    known_landmark_estimates_cds.tags = [estimate["id"],index]
    
    particle_poses = np.array(robot_estimation_trace["timesteps"][index]["particle_poses"])
    update_particle_dict = {
    "x": particle_poses[:,0],
    "y": particle_poses[:,1],
    "weight": estimate["weights"],
    "radius": np.array(estimate["weights"])*50.0,
    }
    particle_poses_cds.data = update_particle_dict
    
    ground_truth_step = ground_truth["timesteps"][ids[index]]
    
    robot_pose_dict = {
    "x": [estimate["robot_pose"][0]],
    "y": [estimate["robot_pose"][1]],
    "degree": [180*estimate["robot_pose"][2]/np.pi],
    "start_angle": [estimate["robot_pose"][2]-opening_angle],
    "end_angle": [estimate["robot_pose"][2]+opening_angle],
    "dist_x_ground_truth": [estimate["robot_pose"][0] - ground_truth_step["robot_pose"][0]],
    "dist_y_ground_truth": [estimate["robot_pose"][1] - ground_truth_step["robot_pose"][1]]}
    
    robot_pose_estimate_cds.data = robot_pose_dict
    landmark_visible = estimate["landmark_visible"]
    landmark_observation_line_dict = {
        "xs" : [[estimate["robot_pose"][0], landmark_pose_trace[index][idf.index(i),0]] for i in landmark_visible], 
        "ys" : [[estimate["robot_pose"][1], landmark_pose_trace[index][idf.index(i),1]] for i in landmark_visible], 
    }
    
    
    landmark_observation_line_cds.data = landmark_observation_line_dict
    
    
    if trace:
        robot_trace_cds.stream({
            "x": [estimate["robot_pose"][0]],
            "y": [estimate["robot_pose"][1]],
            
        })


def update_ground_truth(index, trace=True):
    ground_truth_step = ground_truth["timesteps"][index]
    update_dict = {
        "x": [ground_truth_step["robot_pose"][0]],
        "y": [ground_truth_step["robot_pose"][1]],
    }
    if trace:
        ground_truth_trace_cds.stream(update_dict)
    
    ground_truth_pose_cds.data=update_dict
    
    simulation_state_cds.data["id"]=["Index: %d" % ground_truth_step["id"]],
    simulation_state_cds.data["time"]=["%f sec" % ground_truth_step["timestamp"]]
    


def update(index, trace=False):
    update_observations(index, trace)
    update_ground_truth(ids[index], trace)
    push_notebook(handle=target)

if __name__=="__main__":

    # Precompute landmark trace
    idf = []
    landmark_pose_trace = []
    landmark_ellipse_trace = []

    latest_landmark_poses=[]
    latest_landmark_ellipses = [] #[width, height, rotation]
    times = []
    ids = []

    known_landmark_colors = []
    for estimation in robot_estimation_trace["timesteps"]:
        landmark_colors = ["#%02x%02x%02x" % (200, 100, 0) for i in idf]
        for i,feature_index in enumerate(estimation["landmark_visible"]):
            if not feature_index in idf:
                #new feature
                idf.append(feature_index)
                landmark_colors.append("#%02x%02x%02x" % (255, 0, 0))
                latest_landmark_poses.append(estimation["landmark_poses"][i])
                
                landmark_cov = np.reshape(estimation["landmark_covariances"][i], (2, 2))
                latest_landmark_ellipses.append(
                    covariance_ellipse(landmark_cov))
            else:
                #known feature
                landmark_colors[idf.index(feature_index)] = "#%02x%02x%02x" % (255, 0, 0)
                latest_landmark_poses[idf.index(feature_index)] = estimation["landmark_poses"][i]
                landmark_cov = np.reshape(estimation["landmark_covariances"][i], (2, 2))
                latest_landmark_ellipses[idf.index(feature_index)] =                covariance_ellipse(landmark_cov)
                
        landmark_ellipse_trace.append(np.array(latest_landmark_ellipses.copy()))            
        landmark_pose_trace.append(np.array(latest_landmark_poses.copy()))
        times.append(estimation["timestamp"])
        ids.append(estimation["id"])
        known_landmark_colors.append(landmark_colors)

    #RESET
    landmarks = np.array(ground_truth["landmarks"]) #N_f x 2
    waypoints = np.array(ground_truth["waypoints"]) #N_w x 2

    landmark_colors = ["#%02x%02x%02x" % (0, 0, 255) for i in range(landmarks.shape[0])]
    waypoint_colors = ["#%02x%02x%02x" % (0, 255, 0) for i in range(waypoints.shape[0])]

    landmarks_cds = ColumnDataSource(data = {
        "x": landmarks[:,0],
        "y": landmarks[:,1],
        "color": landmark_colors,
        "landmark_id": range(landmarks.shape[0])
    })

    waypoints_cds = ColumnDataSource(data = {
        "x": waypoints[:,0],
        "y": waypoints[:,1],
        "color": waypoint_colors,
        "waypoint_id": range(waypoints.shape[0])
    })

    ground_truth_step = ground_truth["timesteps"][0]
    ground_truth_trace_cds = ColumnDataSource(data = {
        "x": [ground_truth_step["robot_pose"][0]],
        "y": [ground_truth_step["robot_pose"][1]],
    })

    ground_truth_pose_cds = ColumnDataSource(data = {
        "x": [ground_truth_step["robot_pose"][0]],
        "y": [ground_truth_step["robot_pose"][1]],
    })

    estimate = robot_estimation_trace["timesteps"][0]

    known_landmark_estimates_cds = ColumnDataSource(data = {
        "x": landmark_pose_trace[0],
        "y": landmark_pose_trace[0],
        "color" : known_landmark_colors[0],
        "width": landmark_ellipse_trace[0],
        "height": landmark_ellipse_trace[0],
        "angle": landmark_ellipse_trace[0],
        "degree": landmark_ellipse_trace[0],
        "idf": idf[:len(landmark_pose_trace[0])],
    }, tags = [robot_estimation_trace["timesteps"][0]["id"],0])


    particle_poses = np.array(estimate["particle_poses"])
    particle_poses_cds = ColumnDataSource(data = {
        "x": particle_poses[:,0],
        "y": particle_poses[:,1],
        "weight": estimate["weights"],
        "radius": np.array(estimate["weights"])*50.0,
    })

    opening_angle = np.pi/6
    robot_pose_estimate_cds = ColumnDataSource(data={
        "x": [estimate["robot_pose"][0]],
        "y": [estimate["robot_pose"][1]],
        "degree": [180*estimate["robot_pose"][2]/np.pi],
        "start_angle": [estimate["robot_pose"][2]-opening_angle],
        "end_angle": [estimate["robot_pose"][2]+opening_angle],
        "dist_x_ground_truth": [estimate["robot_pose"][0] - ground_truth_step["robot_pose"][0]],
        "dist_y_ground_truth": [estimate["robot_pose"][1] - ground_truth_step["robot_pose"][1]]
    })

    robot_trace_cds = ColumnDataSource(data = {
        "x": [estimate["robot_pose"][0]],
        "y": [estimate["robot_pose"][1]],
    })

    landmark_visible = estimate["landmark_visible"]
    landmark_observation_line_cds = ColumnDataSource(data = {
        "xs" : [[estimate["robot_pose"][0], landmarks[i,0]] for i in landmark_visible], 
        "ys" : [[estimate["robot_pose"][1], landmarks[i,1]] for i in landmark_visible], 
    })

    simulation_state_cds=ColumnDataSource(data = {
        "x":[-145],
        "y":[90],
        "xid":[-145],
        "yid":[80],
        "id":["Index: %d" % ground_truth_step["id"]],
        "time": ["%f sec" % ground_truth_step["timestamp"]]
    })

    #Plot initialization
    TOOLS="crosshair,pan,wheel_zoom,box_zoom,reset,tap,box_select,lasso_select"

    p = figure(tools=TOOLS, x_range=(-150, 100), y_range=(-100, 100))


    p.axis.major_label_text_font_size = "14pt"
    true_landmarks_artist = p.circle("x","y", radius=1.5, 
                fill_color="color", fill_alpha=0.6, line_color=None, 
                hover_fill_color="black", hover_fill_alpha=0.7, hover_line_color=None, source = landmarks_cds)

    waypoints_line=p.line("x","y",
                line_color="green", line_alpha=0.6, source = waypoints_cds)

    waypoints_artist = p.cross("x","y", size=10, 
                line_color="color", line_alpha=0.6, 
                hover_line_color="black", hover_line_alpha=0.7, source = waypoints_cds)

    landmark_estimate_artist = p.ellipse("x","y", fill_color="color", width="width", height="height", angle="angle",
                                        fill_alpha=1.0, line_color=None, 
                hover_fill_color="black", hover_fill_alpha=0.7, hover_line_color=None, source = known_landmark_estimates_cds)

    ground_truth_trace_artist = waypoints_line=p.line("x","y",
                line_color="blue", line_alpha=0.4, source = ground_truth_trace_cds)

    ground_truth_pose_artist = p.circle("x","y", radius=2, 
                fill_color="cyan", fill_alpha=0.6, line_color=None, 
                hover_fill_color="black", hover_fill_alpha=0.7, hover_line_color=None, source = ground_truth_pose_cds)


    # simulation_time_artist = p.text("x","y", text="time", source = simulation_state_cds)

    # simulation_index_artist = p.text("xid","yid", text="id", source = simulation_state_cds)

    robot_pose_estimate_artist = p.wedge("x","y", radius=12, start_angle="start_angle", end_angle="end_angle",
                fill_color="red", fill_alpha=0.3, line_color=None, 
                hover_fill_color="black", hover_fill_alpha=0.7, hover_line_color=None, source = robot_pose_estimate_cds)

    robot_trace_artist = waypoints_line=p.line("x","y",
                line_color="red", line_alpha=0.4, source = robot_trace_cds)

    landmark_observation_artist = p.multi_line("xs","ys", line_color="red", source = landmark_observation_line_cds)

    particle_pose_artist = p.circle("x","y", radius="radius", 
                fill_color="black", fill_alpha=0.6, line_color=None, 
                hover_fill_color="black", hover_fill_alpha=0.7, hover_line_color=None, source = particle_poses_cds)

    # landmark_hover = HoverTool(renderers = [true_landmarks_artist],tooltips=
    #                         [("Landmark","(@x, @y)"), ("id","@{landmark_id}")])

    # waypoint_hover = HoverTool(renderers = [waypoints_artist],tooltips=
    #                         [("Waypoint","(@x, @y)"), ("id","@{waypoint_id}")])

    # landmark_estimate_hover = HoverTool(renderers = [landmark_estimate_artist],tooltips=
    #                         [("Landmark","(@x, @y)"),
    #                             ("Covariance (width, height,angle)","(@width, @height,@degree)"),
    #                             ("id","@{idf}")])

    # robot_pose_estimate_hover = HoverTool(renderers = [robot_pose_estimate_artist],tooltips=
    #                         [("Position","(@x, @y)"),
    #                             ("Angle","@degree"),
    #                             ("Offset ground truth","(@dist_x_ground_truth,@dist_y_ground_truth)")])

    # p.add_tools(landmark_hover)
    # p.add_tools(waypoint_hover)
    # p.add_tools(landmark_estimate_hover)
    # p.add_tools(robot_pose_estimate_hover)

    driver = webdriver.PhantomJS(executable_path = """C:/Users/phill/node_modules/phantomjs-prebuilt/lib/phantom/bin/phantomjs.exe""")

    # export_png(p, filename="testvis.png", webdriver=driver)
    # image = export_png(p, filename= "testvis.png",webdriver = driver)

    # get and explicit handle to update the next show cell with
    #target = show(p)


    # interact(update, index=widgets.IntSlider(min=0, max=len(ids)-1, step=1, value=0));


    # Run simulation
    video_name = "video.mp4"

    p.toolbar.logo = None
    p.toolbar_location = None
    
    frames=[]
    print("Creating frames...")
    for observation_index, index in enumerate(ids):
        if observation_index==0:
            continue
        if observation_index%10 !=0:
            continue
        print(observation_index)
        update_observations(observation_index)
        update_ground_truth(index)
        frames.append(get_screenshot_as_png(p, driver=driver, width=800, height=800))
        # push updates to the plot continuously using the handle (intererrupt the notebook kernel to stop)
        #push_notebook(handle=target)
    print("Found", len(frames), "frames...")
    print("Saving video as", video_name, "...")
    save_video(frames, video_name)
    

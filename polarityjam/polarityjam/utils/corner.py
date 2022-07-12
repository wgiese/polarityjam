import math

import cv2
import numpy as np
from scipy.spatial import ConvexHull


def get_corner(img, epsilon=5):
    contours, _ = cv2.findContours(img.astype(bool).astype(np.uint8), cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    if len(contours) > 1:
        raise RuntimeError("Too many contours found. Expected a single contour for image!")
    contours = contours[0].squeeze(1)

    hull = ConvexHull(contours)

    corners = douglas_peucker(hull.points, epsilon)

    return corners


def douglas_peucker(points, epsilon):
    if len(points) < 3:
        return points

    dmax = 0
    index = -1
    end = len(points) - 1

    for idx in range(1, end):
        d = perpendicular_distance(points[idx], points[0], points[end])

        if d > dmax:
            index = idx
            dmax = d

    if dmax > epsilon:
        rec_result1 = douglas_peucker(points[0:index + 1], epsilon)
        rec_result2 = douglas_peucker(points[index:], epsilon)

        result = np.concatenate((rec_result1[0:len(rec_result1) - 1], rec_result2), axis=0)

        return result

    else:
        return np.array([points[0], points[end]])


def perpendicular_distance(p, p1, p2):
    if p1[0] == p2[0]:
        result = abs(p[0] - p1[0])
    else:
        slope = (p2[1] - p1[1]) / (p2[0] - p1[0])
        intercept = p1[1] - (slope * p1[0])

        result = abs(slope * p[0] - p[1] + intercept) / math.sqrt(math.pow(slope, 2) + 1)

    return result

from shapely import affinity
from shapely.geometry import Point, box, Polygon, LineString, MultiLineString, LinearRing
from shapely.geometry.base import BaseGeometry

from settings import DASH_SEARCH_BOX_W, DASH_SEARCH_BOX_H, MAX_DOT_LEN, MAX_SWAMP_SYMBOL_LEN, \
    MAX_DASH_LINE_LEN, IMAGE_BBOX_BUFFER, ENDPOINT_FILTER_R, VDOT_FILTER_R, MAX_D2D_DISTANCE, SEARCH_STEP_DEGREE, \
    SML_DASH_SEARCH_BOX_W, SML_DASH_SEARCH_BOX_H, MIN_SWAMP_ENDPOINTS
from line_proc import get_line_endpoints, create_centerline, \
    get_extrapolated_point, get_shortest_path_line, find_nearest_geom, plot_line, \
    get_common_endpoints, get_close_points, filter_geoms, get_search_box, get_points_on_line
from shapely.ops import unary_union, linemerge, nearest_points
from itertools import combinations
from shapely.strtree import STRtree
import matplotlib.pyplot as plt
import logging
from traceback import print_exc
from typing import Tuple, List, Dict


class Dash:
    all_dashes = {}  # key: id(Polygon), value: Dash instance

    # check whether the polygon is already created as Dash
    @staticmethod
    def is_already_created(poly: Polygon) -> bool:
        """
        Check whether this polygon is already created as Dash object.

        :param poly: a Shapely Polygon object
        :return: a boolean value representing whether the polygon is already created as Dash object
        """
        if not isinstance(poly, Polygon):
            raise TypeError(f'Inappropriate type: {type(poly)} for poly whereas a Polygon is expected')
        return id(poly) in Dash.all_dashes

    @staticmethod
    def get_dash_obj(poly: Polygon) -> object:
        """
        Get the Dash object created from the input polygon before.

        :param poly: a Shapely Polygon object
        :return: a Dash object
        """
        if not isinstance(poly, Polygon):
            raise TypeError(f'Inappropriate type: {type(poly)} for poly whereas a Polygon is expected')
        return Dash.all_dashes[id(poly)]

    def __init__(self, poly: Polygon, centerline: BaseGeometry, endpoints: List[Point], dash_line: BaseGeometry = None):
        """
        :param poly: a Shapely Polygon object corresponding to this Dash object
        :param centerline: a Shapely LineString or MultiLineString object created from the polygon corresponding to this
            Dash object
        :param endpoints: a list of Shapely Point objects representing endpoints of the line representing this Dash object
        :param dash_line: a Shapely LineString or MultiLineString object representing the final line representing
            this dash object. The final line is the line drawn for dash itself and its connecting line to its nearby dots.
        """
        if not isinstance(poly, Polygon):
            raise TypeError(f'Inappropriate type: {type(poly)} for poly whereas a Polygon is expected')

        if id(poly) in Dash.all_dashes:
            raise Exception(f'the polygon {poly} already created as Dash object before')

        Dash.all_dashes[id(poly)] = self
        self.poly = poly
        self.centerline = centerline
        self.endpoints = endpoints
        self.dash_line = dash_line
        self.dot_ep_pairs = []

    def save_dot_ep_pairs(self, dot_ep_pairs: List[Tuple[object, Point]]) -> None:
        """
        Store a list of pairs of a Dot object and its closest endpoint of this Dash's centerline.
        These pairs will be used when the connecting lines between this dash and nearby dots are drawn.

        :param dot_ep_pairs: a list of pairs of a Dot object and a Shapely Point object
        :return: None
        """
        for dot, ep in dot_ep_pairs:
            if not isinstance(dot, Dot):
                raise TypeError(f'Inappropriate type: {type(dot)} for dot whereas a Dot is expected')
            if dot in self.dot_ep_pairs:
                raise Exception("The dot is already in this dash")
            if not isinstance(ep, Point):
                raise TypeError(f'Inappropriate type: {type(ep)} for ep whereas a Point is expected')

            self.dot_ep_pairs.append((dot, ep))


class Dot:
    all_dots = {}  # key: id(Point), value: Dot instance

    # check whether the point is already created as Dot
    @staticmethod
    def is_already_created(point: Point) -> bool:
        """
        Check whether this point is already created as Dot object.

        :param point: a Shapely Point object
        :return: a boolean value representing whether the point is already created as Dash object
        """

        if not isinstance(point, Point):
            raise TypeError(f'Inappropriate type: {type(point)} for point whereas a Point is expected')
        return id(point) in Dot.all_dots

    @staticmethod
    def get_dot_obj(point: Point) -> object:
        """
        Get the Dot object created from the input point before.

        :param point: a Shapely Point object
        :return: a Dot object
        """

        if not isinstance(point, Point):
            raise TypeError(f'Inappropriate type: {type(point)} for point whereas a Point is expected')
        return Dot.all_dots[id(point)]

    def __init__(self, point: Point):
        """
        :param point: a Shapely Point object corresponding to this Dot object
        """

        if not isinstance(point, Point):
            raise TypeError(f'Inappropriate type: {type(point)} for point whereas a Point is expected')

        if id(point) in Dot.all_dots:
            raise Exception(f'the point {point} already created as Dot object before')

        Dot.all_dots[id(point)] = self
        self.point = point
        self.dashes = {}

    def search_dashes(self, non_dot_polys_tree: STRtree, poly_line_dict: Dict[int, BaseGeometry], line_ep_dict: Dict[int, List],
                      step_degree: int = SEARCH_STEP_DEGREE, sbox_w: float = DASH_SEARCH_BOX_W, sbox_h: float = DASH_SEARCH_BOX_H)\
            -> List[object]:
        """
        Search and saves the dashes on both sides of the dot
        Returns the searched dashes as tuples

        :param non_dot_polys_tree: STRTree that contains all polygons except the polygons created by dots
        :param poly_line_dict: Dictionary that contains polygon and its created line. {id(polygon):line}.
            This is to prevent redundant calls of creating centerline of polygon.
        :param line_ep_dict: Dictionary that contains line and its endpoints. {id(line):[endpoints]}.
            This is to prevent redundant calls of get_line_endpoints.
        :param step_degree: the rotation degree that will be applied to search boxes after each search iteration
        :param sbox_h: height of dash search box
        :param sbox_w: width of dash search box
        :return: a list of tuples that contains attributes for a dash (poly, centerline, target_ep, endpoints)
        """

        def get_dash_search_boxes(dot: object, width: float = sbox_w, height: float = sbox_h) -> List[object]:
            """
            Get two dash search boxes. Each box covers one side from the dot.
            (e.g., box1: left side of dot, box2: right side of dot)
            The reason it creates two separate boxes is to reduce detecting false dashes.
            The search algorithm regards polygons as dashes only if both boxes search exactly one polygon each.
            The size of these boxes are assumed to be the distance between this dot and nearby dashes.

            :param dot: a Dot object that will be the center of the search boxes
            :param width: width of dash search box. This value corresponds to the widths of both boxes.
            :param height: height of dash search box. This value corresponds to the heights of both boxes.
            :return: a list of Shapely Polygon objects
            """
            box1 = box(dot.x - width / 2, dot.y - height / 2, dot.x, dot.y + height / 2)
            box2 = box(dot.x, dot.y - height / 2, dot.x + width / 2, dot.y + height / 2)
            return [box1, box2]

        # the dash polygon and endpoint pairs found in the entire search process of a dot
        found_poly_ep_pairs = set()  # elements: f'{str(id(Polygon))}{str(id(Point))}'
        searched_dashes_all = []  # all the dash tuples created by the entire search process of a dot

        default_boxes = get_dash_search_boxes(self.point)

        # rotate the search boxes to find the nearby two dashes
        for d in range(0, 180, step_degree):
            search_boxes = [affinity.rotate(sbox, d, origin=self.point) for sbox in default_boxes]

            all_found = True  # True when all search boxes found exactly one polygon each
            searched_dashes_two_boxes = []  # the polygons searched by this pair of search boxes

            # handle each search box one at a time
            for sbox in search_boxes:
                if not all_found:
                    continue

                # candidate polygons searched by one search box
                candidate_polys = [geom for geom in non_dot_polys_tree.query(sbox)
                                   if geom.intersects(sbox)]

                # filter out false dash polygons based on rules
                # dashes searched in a box; [(polygon, centerline, endpoints)]
                searched_dashes_one_box = []  # searched dashes in one box
                for poly in candidate_polys:
                    # rule 1: dash polygon's perimeter should long enough
                    if poly.length > MAX_DOT_LEN:
                        # get centerline of the polygon
                        poly_centerline = create_centerline(poly, poly_line_dict)

                        # get endpoints of the centerline
                        endpoints = get_line_endpoints(poly_centerline, line_ep_dict)

                        # rule 2: dash polygon is not swamp symbol (too many endpoints compared to its perimeter)
                        if not (len(endpoints) > MIN_SWAMP_ENDPOINTS and poly.length < MAX_SWAMP_SYMBOL_LEN):
                            # rule 3: dash polygon's endpoint should be in the search box
                            tree = STRtree(endpoints)
                            endpoints_in_sbox = tree.query(sbox)
                            endpoints_in_sbox = [ep for ep in endpoints_in_sbox if sbox.covers(ep)]
                            if len(endpoints_in_sbox) > 0:
                                # find the nearest endpoint
                                target_ep = find_nearest_geom(self.point, endpoints_in_sbox)

                                # filter out false endpoints that is assumed to be created
                                # at branches near the actual endpoint
                                filtered_endpoints = filter_geoms(target_ep, endpoints, ENDPOINT_FILTER_R)
                                filtered_endpoints.append(target_ep)

                                temp_dash = (poly, poly_centerline, target_ep, filtered_endpoints)
                                searched_dashes_one_box.append(temp_dash)

                # rule 4: one search box should find only one dash
                # however, the dashes found by each box can be the same dash
                # (in case there are other lines or symbols overlaps with two nearby dashes at the same time)
                if len(searched_dashes_one_box) == 1:
                    searched_dashes_two_boxes.extend(searched_dashes_one_box)
                else:
                    all_found = False
                    searched_dashes_two_boxes.clear()

            # for sbox in search_boxes:
            #     plt.plot(*sbox.exterior.xy)
            if all_found:
                # for sbox in search_boxes:
                #     plt.plot(*sbox.exterior.xy)

                for dash_tuple in searched_dashes_two_boxes:
                    poly, _, target_ep, _ = dash_tuple
                    poly_ep_id = f'{str(id(poly))}{str(id(target_ep))}'

                    # add to searched dashes list only if it is a new polygon-endpoint pair
                    if poly_ep_id not in found_poly_ep_pairs:
                        found_poly_ep_pairs.add(poly_ep_id)

                        searched_dashes_all.append(dash_tuple)

        return searched_dashes_all

    def search_additional_dashes(self, non_dot_polys_tree: STRtree, step_degree: int = SEARCH_STEP_DEGREE,
                                 sbox_w: float = SML_DASH_SEARCH_BOX_W, sbox_h: float = SML_DASH_SEARCH_BOX_H) -> List[object]:
        """
        Search and saves additional dashes around the dot.
        Only the polygons that are already detected as dashes by other dots are searched.
        Returns the additionally searched and updated Dash objects

        :param non_dot_polys_tree: STRTree that contains all polygons except the polygons created by dots
        :param step_degree: the rotation degree that will be applied to search boxes after each search iteration
        :param sbox_h: height of dash search box
        :param sbox_w: width of dash search box
        :return: a list of additionally searched dash polygons
        """

        default_box = get_search_box(self.point, sbox_w, sbox_h)
        searched_polys = []
        searched_poly_ids = set()

        # rotate the search box to find the nearby ADDITIONAL dashes
        for d in range(0, 180, step_degree):
            sbox = affinity.rotate(default_box, d, origin=self.point)
            # plt.plot(*sbox.exterior.xy)

            # candidate polygons searched by the search box
            # filter out polygons that are already searched by previous iteration of this search process
            # filter out polygons that are already associated with this dot
            # filter out polygons that are not created as a Dash object before
            for geom in non_dot_polys_tree.query(sbox):
                if (geom.intersects(sbox)
                        and id(geom) not in searched_poly_ids
                        and id(geom) not in self.dashes
                        and id(geom) in Dash.all_dashes):
                    searched_polys.append(geom)
                    searched_poly_ids.add(id(geom))

        return searched_polys

    def associate_dashes(self, non_dot_polys_tree: STRtree, poly_line_dict: Dict[int, BaseGeometry], line_ep_dict: Dict[int, List]) \
            -> List[object]:
        """
        Search dashes and associate searched dashes with the current dot object.
        The searched dash objects are created. The current dot object is also associated with the created dash objects

        :param non_dot_polys_tree: STRTree that contains all polygons except the polygons created by dots
        :param poly_line_dict: Dictionary that contains polygon and its created line. {id(polygon):line}.
            This is to prevent redundant calls of creating centerline of polygon.
        :param line_ep_dict: Dictionary that contains line and its endpoints. {id(line):[endpoints]}.
            This is to prevent redundant calls of get_line_endpoints.
        :return: a list of dash objects
        """

        dash_tuples = self.search_dashes(non_dot_polys_tree, poly_line_dict, line_ep_dict)

        associated_dashes = []
        # create or update dash objects based on searched dash tuples around the dot
        for poly, centerline, target_ep, endpoints in dash_tuples:
            # update a Dash object
            if Dash.is_already_created(poly):
                dash = Dash.get_dash_obj(poly)
                dash.save_dot_ep_pairs([(self, target_ep)])
                # if dash object has already made by other dots, then store only the common endpoints
                dash.endpoints = get_common_endpoints(dash.endpoints, endpoints)
            # create a new Dash object
            else:
                dash = Dash(poly, centerline, endpoints)
                dash.save_dot_ep_pairs([(self, target_ep)])

            self.dashes[id(poly)] = dash
            associated_dashes.append(dash)

        return associated_dashes


def extract_dot_dashed_lines(dot_ps: List[Point], polygons: List[Polygon], image_bbox: LinearRing) -> BaseGeometry:
    """
    Extracts a solid vector line on dot-dash lines on the map.
    While searching dots and dashes on the map, it also creates
    instances of Dot and Dash classes.

    :param dot_ps: dots as shapely Points
    :param polygons: all the polygons on the map as shapely Polygons
    :param image_bbox: image bounding box as a Shapely LinearRing object
    :return: connected dot-dashed lines as a Shapely MultiString
    """

    all_drawn_lines = []
    poly_line_dict = {}  # id(polygon):centerline, to prevent creating centerline for the same polygons
    line_ep_dict = {}
    dot_ps_tree = STRtree(dot_ps)

    # filter out polygons created by dots; the perimeters of polygons are used to distinguish dot polygons
    dot_ps_to_remove = set()  # {id(Point), }
    non_dot_polys = [poly for poly in polygons if poly.length > MAX_DOT_LEN]

    # for poly in non_dot_polys:
    #     plt.plot(*poly.exterior.xy)

    # filter out dot points that are within non dot polygons (considered as false dot points)
    for poly in non_dot_polys:
        overlapped_dot_ps = dot_ps_tree.query(poly)
        dot_ps_to_remove.update([id(dot_p) for dot_p in overlapped_dot_ps if poly.covers(dot_p)])
    dot_ps = [dot_p for dot_p in dot_ps if id(dot_p) not in dot_ps_to_remove]

    # for dot_p in dot_ps:
    #     plt.plot(dot_p.x, dot_p.y, marker="D")

    non_dot_polys_tree = STRtree(non_dot_polys)

    logging.info('Searching dash polygons around the detected dots')
    # initial dot and dash object creation
    for p in dot_ps:
        dot = Dot(p)
        dot.associate_dashes(non_dot_polys_tree, poly_line_dict, line_ep_dict)

    logging.info('Making virtual dots and searching dash polygons around them')
    # make virtual dots (dots not detected) for dashes having dots less than two.
    dashes_copy = list(Dash.all_dashes.values())
    vdot_ps = []  # a list of Shapely Point objects representing virtual dots' locations
    for dash in dashes_copy:
        if dash.centerline is not None:
            # plot_line(dash.centerline)
            # filter out the endpoints close to valid dots
            # because virtual dot is no need at that end
            endpoints_wo_dots = dash.endpoints.copy()

            for _, ep in dash.dot_ep_pairs:
                if ep in endpoints_wo_dots:
                    endpoints_wo_dots.remove(ep)

            # get the locations of virtual dots using extrapolation
            for ep in endpoints_wo_dots:
                vdot_p = get_extrapolated_point(dash.centerline, ep)
                vdot_ps.append(vdot_p)
                # plt.plot(vdot_p.x, vdot_p.y, marker="*")

    # make virtual dots
    redundant_vdot_ps = set()  # {id(Point),}
    vdots_tree = STRtree(vdot_ps)

    # find out redundant virtual dot points that are too close to the valid dots
    for dot in Dot.all_dots.values():
        filter_circle = dot.point.buffer(VDOT_FILTER_R)
        filtered = vdots_tree.query(filter_circle)
        redundant_vdot_ps.update([id(p) for p in filtered if filter_circle.covers(p)])
        # plt.plot(dot.point.x, dot.point.y, marker="o")
        # plt.plot(*filter_circle.exterior.xy)

    # create dot objects of virtual dot points
    for vdot_p in vdot_ps:
        # create only the virtual dots that are not redundant
        if id(vdot_p) not in redundant_vdot_ps:
            vdot = Dot(vdot_p)
            associated_dashes = vdot.associate_dashes(non_dot_polys_tree, poly_line_dict, line_ep_dict)

            # remove virtual dot when there is no dash around it (false virtual dot)
            if len(associated_dashes) == 0:
                # plt.plot(vdot_p.x, vdot_p.y, marker="x")
                del Dot.all_dots[id(vdot.point)]
                del vdot  # delete virtual dot object
            else:
                # filter out redundant virtual dot points that are too close to the dot just created
                # plt.plot(vdot_p.x, vdot_p.y, marker="*")
                filter_circle = vdot_p.buffer(VDOT_FILTER_R)
                filtered = vdots_tree.query(filter_circle)
                redundant_vdot_ps.update([id(p) for p in filtered if filter_circle.covers(p)])

    # search and associate additional dashes
    for dot in Dot.all_dots.values():
        dash_polys = dot.search_additional_dashes(non_dot_polys_tree)

        # update the dash objects
        # associate the searched dash with this dot
        for poly in dash_polys:
            dash = Dash.get_dash_obj(poly)
            # plot_line(dash.centerline)
            all_points_on_line = get_points_on_line(dash.centerline)

            contact_point = nearest_points(all_points_on_line, dot.point)[0]  # contact point on the centerline

            # update the dash
            dash.save_dot_ep_pairs([(dot, contact_point)])

            # store the updated dash in the dot object
            dot.dashes[id(poly)] = dash

    # obtain dash lines
    logging.info('Extracting the lines from dots and dash polygons')
    buffered_image_bbox = image_bbox.buffer(IMAGE_BBOX_BUFFER, single_sided=True)
    for dash in Dash.all_dashes.values():
        if dash.centerline is not None:
            # plot_line(dash.centerline)
            # prepare dash body line to be merged
            if isinstance(dash.centerline, MultiLineString):
                lines = list(dash.centerline.geoms)
            else:
                lines = [dash.centerline]

            # make the connecting line from dash's endpoint to dot
            for dot, ep in dash.dot_ep_pairs:
                lines.append(LineString([ep, dot.point]))

            # make virtual Dot objects to draw the lines connecting between dash to image bounding box
            # points that intersect with buffered image bounding box
            close_points = get_close_points(dash.centerline, buffered_image_bbox.interiors[0])
            for cp in close_points:
                p_bbox = nearest_points(image_bbox, cp)[0]
                # create a virtual dot on the image bounding box and register it to dash
                dash.dot_ep_pairs.append((Dot(p_bbox), cp))
                lines.append(LineString([cp, p_bbox]))

            # merged line with the dash body line and connecting line to dots
            mline = linemerge(lines)
            # to split lines at every intersection; this is to create intersection vertices in graph
            mline = unary_union(mline)

            # find the shortest paths between dots
            logging.info('Finding shortest path between dots to connect dots and dashes')
            path_lines = []
            for pair1, pair2 in combinations(dash.dot_ep_pairs, 2):
                # if the points are too far away from each other, it is likely a solid line in between
                if pair1[0].point.distance(pair2[0].point) < MAX_D2D_DISTANCE:
                    path_line = get_shortest_path_line(mline, pair1[0].point, pair2[0].point)
                    path_lines.append(path_line)

            # filter out drawn lines that are possibly solid lines
            logging.info('Filtering out long lines')
            path_lines = [line for line in path_lines
                          if line is not None and line.length < MAX_DASH_LINE_LEN]

            if len(path_lines) > 0:
                logging.info('Adding extracted dot dash lines')
                try:
                    dash_line = unary_union(path_lines)
                    dash.dash_line = dash_line

                    if isinstance(dash_line, MultiLineString):
                        all_drawn_lines.extend(list(dash_line.geoms))
                    elif isinstance(dash_line, LineString):  # LineString
                        all_drawn_lines.append(dash_line)
                    else:
                        TypeError(
                            f'Inappropriate type: {type(dash_line)} for dash_line whereas a MultiLineString or LineString is expected')

                except Exception as err:
                    logging.error(err, exc_info=True)
                    print_exc()
                    print(err)

    logging.info('Adding image bounding lines')
    all_drawn_lines.append(image_bbox)
    all_drawn_lines = unary_union(all_drawn_lines)

    logging.info('Merging all extracted lines')
    return linemerge(all_drawn_lines)
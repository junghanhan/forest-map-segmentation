from shapely import affinity
from shapely.geometry import Point, box, Polygon, LineString, MultiLineString, MultiPoint
from settings import DASH_SEARCH_BOX_W, DASH_SEARCH_BOX_H, MAX_DOT_LEN, MAX_SWAMP_SYMBOL_LEN, \
    MAX_DASH_LINE_LEN, IMAGE_BBOX_BUFFER, ENDPOINT_FILTER_R, VDOT_FILTER_R, MAX_P2P_DISTANCE, SEARCH_STEP_DEGREE
from line_proc import get_line_endpoints, create_centerline, \
    get_extrapolated_point, get_path_line, find_nearest_geom, plot_line, \
    get_common_endpoints, get_close_points, filter_geoms, get_search_box, get_points_on_line
from shapely.ops import unary_union, linemerge, nearest_points
from itertools import combinations
from shapely.strtree import STRtree
import matplotlib.pyplot as plt
import logging
from traceback import print_exc


# TODO: getter, setter implement (property)
class Dash:
    all_dashes = {}  # key: id(Polygon), value: Dash instance

    # check whether the polygon is already created as Dash
    @staticmethod
    def is_already_created(poly):
        if not isinstance(poly, Polygon):
            raise TypeError(f'Inappropriate type: {type(poly)} for poly whereas a Polygon is expected')
        return id(poly) in Dash.all_dashes

    @staticmethod
    def get_dash_obj(poly):
        if not isinstance(poly, Polygon):
            raise TypeError(f'Inappropriate type: {type(poly)} for poly whereas a Polygon is expected')
        return Dash.all_dashes[id(poly)]

    def __init__(self, poly, centerline, endpoints, dash_line=None):
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

    def save_dot_ep_pairs(self, dot_ep_pairs):
        for dot, ep in dot_ep_pairs:
            if not isinstance(dot, Dot):
                raise TypeError(f'Inappropriate type: {type(dot)} for dot whereas a Dot is expected')
            if dot in self.dot_ep_pairs:
                raise Exception("The dot is already in this dash")
            if not isinstance(ep, Point):
                raise TypeError(f'Inappropriate type: {type(ep)} for ep whereas a Point is expected')

            self.dot_ep_pairs.append((dot, ep))


# TODO: getter, setter implement (property)
class Dot:
    all_dots = {}  # key: id(Point), value: Dot instance

    # check whether the point is already created as Dot
    @staticmethod
    def is_already_created(point):
        if not isinstance(point, Point):
            raise TypeError(f'Inappropriate type: {type(point)} for point whereas a Point is expected')
        return id(point) in Dot.all_dots

    @staticmethod
    def get_dot_obj(point):
        if not isinstance(point, Point):
            raise TypeError(f'Inappropriate type: {type(point)} for point whereas a Point is expected')
        return Dot.all_dots[id(point)]

    def __init__(self, point):
        if not isinstance(point, Point):
            raise TypeError(f'Inappropriate type: {type(point)} for point whereas a Point is expected')

        if id(point) in Dot.all_dots:
            raise Exception(f'the point {point} already created as Dot object before')

        Dot.all_dots[id(point)] = self
        self.point = point
        self.dashes = {}

    def search_dashes(self, non_dot_polys_tree, poly_line_dict, line_ep_dict,
                      step_degree=SEARCH_STEP_DEGREE, sbox_w=DASH_SEARCH_BOX_W, sbox_h=DASH_SEARCH_BOX_H):
        """
        Search and saves the dashes on both sides of the dot
        Returns the searched Dash objects

        side effect: Dash objects can be created or updated.

        :param non_dot_polys_tree: STRTree that contains all polygons except the polygons created by dots
        :param poly_line_dict: Dictionary that contains polygon and its created line. {id(polygon):line}.
            This is to prevent redundant calls of creating centerline of polygon.
        :param line_ep_dict: Dictionary that contains line and its endpoints. {id(line):[endpoints]}.
            This is to prevent redundant calls of get_line_endpoints.
        :param step_degree: the rotation degree that will be applied to search boxes after each search iteration
        :param sbox_h: height of dash search box
        :param sbox_w: width of dash search box
        :return: a list of dash objects
        """

        # distance: assumed distance from dot to dash (2*distance = search box width)
        # buffer: search buffer (2*buffer = search box height)
        def get_dash_search_boxes(dot, width=sbox_w, height=sbox_h):
            box1 = box(dot.x - width / 2, dot.y - height / 2, dot.x, dot.y + height / 2)
            box2 = box(dot.x, dot.y - height / 2, dot.x + width / 2, dot.y + height / 2)
            return [box1, box2]

        # the dash polygon and endpoint pairs found in the search process
        found_poly_ep_pairs = set()  # elements: f'{str(id(Polygon))}{str(id(Point))}'
        found_dashes = []  # the Dash objects created by the search process

        default_boxes = get_dash_search_boxes(self.point)

        # rotate the search boxes to find the nearby two dashes
        for d in range(0, 180, step_degree):
            search_boxes = [affinity.rotate(sbox, d, origin=self.point) for sbox in default_boxes]

            all_found = True  # True when all search boxes found exactly one polygon each
            all_searched_dashes = []  # the polygons searched by this pair of search boxes

            # handle each search box one at a time
            for sbox in search_boxes:
                if not all_found:
                    continue

                # candidate polygons searched by one search box
                candidate_polys = [geom for geom in non_dot_polys_tree.query(sbox)
                                   if geom.intersects(sbox)]

                # filter out false dash polygons based on rules
                # dashes searched in a box; [(polygon, centerline, endpoints)]
                searched_dashes = []
                for poly in candidate_polys:
                    # rule 1: dash polygon's perimeter should long enough
                    if poly.length > MAX_DOT_LEN:
                        # get centerline of the polygon
                        poly_centerline = create_centerline(poly, poly_line_dict)

                        # get endpoints of the centerline
                        endpoints = get_line_endpoints(poly_centerline, line_ep_dict)

                        # rule 2: dash polygon is not swamp symbol (too many endpoints compared to its perimeter)
                        # TODO: 5 is a temporary value; to filter out swamp symbol with 6 endpoints
                        if not (len(endpoints) > 5 and poly.length < MAX_SWAMP_SYMBOL_LEN):
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
                                searched_dashes.append(temp_dash)

                # rule 4: one search box should find only one dash
                # however, the dashes found by each box can be the same dash
                if len(searched_dashes) == 1:
                    all_searched_dashes.extend(searched_dashes)
                else:
                    all_found = False
                    all_searched_dashes.clear()

            # for sbox in search_boxes:
            #     plt.plot(*sbox.exterior.xy)
            if all_found:
                # for sbox in search_boxes:
                #     plt.plot(*sbox.exterior.xy)

                # create dash objects and pair them
                for poly, centerline, target_ep, endpoints in all_searched_dashes:
                    poly_ep_id = f'{str(id(poly))}{str(id(target_ep))}'

                    # create or update dash object only if it is a new polygon-endpoint pair
                    if poly_ep_id not in found_poly_ep_pairs:
                        found_poly_ep_pairs.add(poly_ep_id)

                        if Dash.is_already_created(poly):
                            dash = Dash.get_dash_obj(poly)
                            dash.save_dot_ep_pairs([(self, target_ep)])
                            # if dash object has already made by other dots, then store only the common endpoints
                            dash.endpoints = get_common_endpoints(dash.endpoints, endpoints)
                        else:
                            dash = Dash(poly, centerline, endpoints)
                            dash.save_dot_ep_pairs([(self, target_ep)])

                        found_dashes.append(dash)

                        # store the created dashes in the dot object
                        self.dashes[id(poly)] = dash

        return found_dashes

    def search_additional_dashes(self, non_dot_polys_tree, step_degree=SEARCH_STEP_DEGREE,
                                 sbox_w=DASH_SEARCH_BOX_W / 1.5, sbox_h=DASH_SEARCH_BOX_H / 1.5):
        """
        Search and saves additional dashes around the dot.
        Only the polygons that are already detected as dashes by other dots are searched.
        Returns the additionally searched and updated Dash objects

        side effect: Dash objects can be updated.

        :param non_dot_polys_tree: STRTree that contains all polygons except the polygons created by dots
        :param step_degree: the rotation degree that will be applied to search boxes after each search iteration
        :param sbox_h: height of dash search box
        :param sbox_w: width of dash search box
        :return: a list of updated dash objects
        """

        default_box = get_search_box(self.point, sbox_w, sbox_h)
        updated_dashes = []

        # rotate the search box to find the nearby ADDITIONAL dashes
        for d in range(0, 180, step_degree):
            sbox = affinity.rotate(default_box, d, origin=self.point)
            # plt.plot(*sbox.exterior.xy)

            # candidate polygons searched by the search box
            # filter out polygons that are already associated with this dot
            # filter out polygons that are not created as a Dash object before
            searched_polys = [geom for geom in non_dot_polys_tree.query(sbox)
                               if geom.intersects(sbox)
                               and id(geom) not in self.dashes
                               and id(geom) in Dash.all_dashes]

            # if len(searched_polys) > 0:
            #     plt.plot(*sbox.exterior.xy)
            # for poly in searched_polys:
            #     plt.plot(*poly.exterior.xy)

            # associate the dash object with this dot
            for poly in searched_polys:
                dash = Dash.get_dash_obj(poly)
                # plot_line(dash.centerline)
                all_points_on_line = get_points_on_line(dash.centerline)

                contact_point = nearest_points(all_points_on_line, self.point)[0]  # contact point on the centerline

                # update the dash
                dash.save_dot_ep_pairs([(self, contact_point)])

                # store the updated dash in the dot object
                self.dashes[id(poly)] = dash

                updated_dashes.append(dash)

        return updated_dashes













# Helper functions related to Dot, Dash object
# deprecated
# def create_dash_pairs(dot, dash_poly_pairs):
#     """
#     Create dash pairs from dash polygon pairs
#     When Dash object is already created from the same polygon, the dot is stored in that Dash object.
#     When Dash object is newly created, the input dot is stored in that Dash object
#
#     :param dot: Dot object
#     :param dash_poly_pairs: dash polygon pairs around the Dot object
#     :return: Dash object pairs around the Dot object
#     """
#
#     dash_pairs = []
#     for dp1, dp2 in dash_poly_pairs:
#         if Dash.is_already_created(dp1):
#             dash1 = Dash.get_dash_obj(dp1)
#             dash1.save_dots([dot])
#         else:
#             dash1 = Dash(dp1, [dot])
#
#         if Dash.is_already_created(dp2):
#             dash2 = Dash.get_dash_obj(dp2)
#             dash2.save_dots([dot])
#         else:
#             dash2 = Dash(dp2, [dot])
#
#         dash_pairs.append((dash1, dash2))
#
#     return dash_pairs


def extract_dot_dashed_lines(dot_ps, polygons, image_bbox):
    """
    Extracts a solid vector line on dot-dash lines on the map.
    While searching dots and dashes on the map, it also creates
    instances of Dot and Dash classes.
    TODO: need refactoring and performance optimization

    :param dot_ps: dots as shapely Points
    :param polygons: all the polygons on the map as shapely Polygons
    :param image_bbox: image bounding box as a Shapely LinearRing object
    :return: connected dot-dashed lines as a Shapely MultiString
    """

    all_drawn_lines = []
    poly_line_dict = {}  # id(polygon):centerline
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
        dashes = dot.search_dashes(non_dot_polys_tree, poly_line_dict, line_ep_dict)

    logging.info('Making virtual dots and searching dash polygons around them')
    # make virtual dots (dots not detected) for dashes having dots less than two.
    dashes_copy = list(Dash.all_dashes.values())
    vdot_ps = []  # a list of Shapely Point objects representing virtual dots' locations
    for dash in dashes_copy:
        # TODO: need to come up with a way to determine if a dash has an endpoint without dot associated
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
            dashes = vdot.search_dashes(non_dot_polys_tree, poly_line_dict, line_ep_dict)
            if len(dashes) == 0:
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
        dot.search_additional_dashes(non_dot_polys_tree)

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
                if pair1[0].point.distance(pair2[0].point) < MAX_P2P_DISTANCE:
                    path_line = get_path_line(mline, pair1[0].point, pair2[0].point)
                    path_lines.append(path_line)

            # filter out drawn lines that are possibly solid lines
            logging.info('Filtering out long lines')
            path_lines = [line for line in path_lines
                          if line is not None and line.length < MAX_DASH_LINE_LEN]

            if len(path_lines) > 0:
                logging.info('Adding extracted dot dash lines')
                # TODO: temporary exception handling for debugging purpose
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
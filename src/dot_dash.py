from shapely import affinity
from shapely.geometry import Point, box, Polygon, LineString, MultiLineString
from settings import DASH_SEARCH_BOX_W, DASH_SEARCH_BOX_H, MAX_DOT_LENGTH, MAX_SWAMP_SYMBOL_LEN
from line_proc import is_endpoint_inside, get_line_endpoints, create_centerline, \
    get_extrapolated_point, get_path_line, find_nearest_geom, plot_line
from shapely.ops import unary_union, linemerge
from itertools import combinations
from shapely.strtree import STRtree
import matplotlib.pyplot as plt
import logging

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

    def __init__(self, poly, conn_dots=[], dash_line=None):
        if not isinstance(poly, Polygon):
            raise TypeError(f'Inappropriate type: {type(poly)} for poly whereas a Polygon is expected')

        if id(poly) in Dash.all_dashes:
            raise Exception(f'the polygon {poly} already created as Dash object before')

        Dash.all_dashes[id(poly)] = self
        self.poly = poly
        self.conn_dots = []
        self.save_dots(conn_dots)
        self.dash_line = dash_line

    def save_dots(self, dots):
        for dot in dots:
            if not isinstance(dot, Dot):
                raise TypeError(f'Inappropriate type: {type(dot)} for dot whereas a Dot is expected')
            if dot in self.conn_dots:
                raise Exception("The dot is already in this dash")
            self.conn_dots.append(dot)

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

    def __init__(self, point, dash_pairs=[]):
        if not isinstance(point, Point):
            raise TypeError(f'Inappropriate type: {type(point)} for point whereas a Point is expected')

        if id(point) in Dot.all_dots:
            raise Exception(f'the point {point} already created as Dot object before')

        Dot.all_dots[id(point)] = self
        self.point = point
        self.dashes = {}
        self.dash_pairs = []
        self.save_dash_pairs(dash_pairs)

    def search_dash_polygons(self, strtree, poly_line_dict, step_degree=20):
        """
        Search the dashes on both sides of the dot
        Returns the pairs of polygons that are assumed as dashes

        strtree: STRTree that contains all polygons except the polygons created by dots
        poly_line_dict: Dictionary that contains polygon and its created line. {id(polygon):line}.
            This is to prevent redundant calls of creating centerline of polygon.
        """

        # distance: assumed distance from dot to dash (2*distance = search box width)
        # buffer: search buffer (2*buffer = search box height)
        def get_dash_search_boxes(dot, width=DASH_SEARCH_BOX_W, height=DASH_SEARCH_BOX_H):
            box1 = box(dot.x - width / 2, dot.y - height / 2, dot.x, dot.y + height / 2)
            box2 = box(dot.x, dot.y - height / 2, dot.x + width / 2, dot.y + height / 2)
            return [box1, box2]

        found_dashes = {}
        dash_pairs = []

        default_boxes = get_dash_search_boxes(self.point)

        # rotate the search boxes to find the nearby two dashes
        for d in range(0, 180, step_degree):
            search_boxes = [affinity.rotate(sbox, d, origin=self.point) for sbox in default_boxes]

            all_found = True  # True when all search boxes found exactly one polygon each
            all_searched_polys = []  # the polygons searched by all the search boxes rotated with a certain degree

            # handle each search box one at a time
            for sbox in search_boxes:
                if not all_found:
                    continue

                # searched polygons by one search box
                candidate_polys = [geom for geom in strtree.query(sbox) if
                                  geom.intersects(sbox) and not id(geom) in found_dashes]

                # check if the found polygon's endpoint is within the search box
                # and its length is long enough to be a dash
                # TODO: How about storing centerline instead of polygons in the first place?
                searched_polys = []
                for poly in candidate_polys:
                    if poly.length > MAX_DOT_LENGTH:
                        if id(poly) in poly_line_dict:
                            poly_centerline = poly_line_dict[id(poly)]
                        else:
                            poly_centerline = create_centerline(poly)
                            poly_line_dict[id(poly)] = poly_centerline

                        # filter out polygons
                        endpoints = get_line_endpoints(poly_centerline)
                        # if at least an endpoint is in the search box
                        if any(sbox.contains(p) for p in endpoints):
                            # if not swamp symbol
                            # TODO: 5 is a temporary value; to filter out swamp symbol with 6 endpoints
                            if not(len(endpoints) > 5 and poly.length < MAX_SWAMP_SYMBOL_LEN):
                                searched_polys.append(poly)

                # check if each search box found only one
                if len(searched_polys) == 1 and not searched_polys[0] in all_searched_polys:
                    all_searched_polys.extend(searched_polys)
                else:
                    all_found = False
                    all_searched_polys.clear()

            if all_found:
                for poly in all_searched_polys:
                    found_dashes[id(poly)] = poly
                if len(all_searched_polys) != 2:
                    raise Exception("all_searched_polys does not contain 2 polygons")
                dash_pairs.append(tuple(all_searched_polys))

                # plot
                # for sbox in search_boxes:
                #     bbox_x, bbox_y = sbox.exterior.xy
                #     plt.plot(bbox_x, bbox_y)

        return dash_pairs

    def save_dash_pairs(self, dash_pairs):
        for dash_pair in dash_pairs:
            for dash in dash_pair:
                if not isinstance(dash, Dash):
                    raise TypeError(f'Inappropriate type: {type(dash)} for dash whereas a Dash is expected')
                if dash in self.dashes:
                    raise Exception("The dash is already in this dot")
                self.dashes[id(dash)] = dash
            self.dash_pairs.append(dash_pair)


# Helper functions related to Dot, Dash object
def create_dash_pairs(dot, dash_poly_pairs):
    """
    Create dash pairs from dash polygon pairs
    When Dash object is already created from the same polygon, the dot is stored in that Dash object.
    When Dash object is newly created, the input dot is stored in that Dash object

    :param dot: Dot object
    :param dash_poly_pairs: dash polygon pairs around the Dot object
    :return: Dash object pairs around the Dot object
    """

    dash_pairs = []
    for dp1, dp2 in dash_poly_pairs:
        if Dash.is_already_created(dp1):
            dash1 = Dash.get_dash_obj(dp1)
            dash1.save_dots([dot])
        else:
            dash1 = Dash(dp1, [dot])

        if Dash.is_already_created(dp2):
            dash2 = Dash.get_dash_obj(dp2)
            dash2.save_dots([dot])
        else:
            dash2 = Dash(dp2, [dot])

        dash_pairs.append((dash1, dash2))

    return dash_pairs


def extract_dot_dashed_lines(dots, polygons, max_dot_length=0.0005):
    """
    Extracts a solid vector line on dot-dash lines on the map.
    While searching dots and dashes on the map, it also creates
    instances of Dot and Dash classes.
    TODO: need refactoring and performance optimization

    :param dots: dots as shapely Points
    :param polygons: all the polygons on the map as shapely Polygons
    :param max_dot_length: maximum length of dot Shapely Polygon; used to ignore dot polygons
        when dash polygons are searched
    :return: connected dot-dashed lines as a Shapely MultiString
    """

    all_drawn_lines = []
    poly_line_dict = {}  # id(polygon):centerline
    dots_tree = STRtree(dots)
    polygons_wo_dots = [poly for poly in polygons if len([dot for dot in dots_tree.query(poly) if
                                                          poly.intersects(dot) and poly.length <= max_dot_length]) == 0]

    # plotting non dot polygons
    # for geom in polygons_wo_dots:
    #   x,y = geom.exterior.xy # x,y are arrays
    #   plt.plot(x,y)

    polygons_wo_dots_tree = STRtree(polygons_wo_dots)

    logging.info('Searching dash polygons around the detected dots')
    for p in dots:
        dot = Dot(p)
        dash_poly_pairs = dot.search_dash_polygons(polygons_wo_dots_tree, poly_line_dict)
        if len(dash_poly_pairs) == 0:
            del Dot.all_dots[id(p)]
        else:
            dash_pairs = create_dash_pairs(dot, dash_poly_pairs)  # Dash object is created
            dot.save_dash_pairs(dash_pairs)

    logging.info('Making virtual dots and searching dash polygons around them')
    # make virtual dots (dots not detected) for dashes having dots less than two
    dashes_copy = list(Dash.all_dashes.values())
    for dash in dashes_copy:
        if len(dash.conn_dots) < 2:
            if id(dash.poly) in poly_line_dict:
                dash_body_line = poly_line_dict[id(dash.poly)]
            else:
                dash_body_line = create_centerline(dash.poly)
                poly_line_dict[id(dash.poly)] = dash_body_line

            if dash_body_line is not None:
                endpoints = get_line_endpoints(dash_body_line)

                # filter out the endpoints close to valid dots
                # because virtual dot is no need at that end
                endpoints_wo_dots = endpoints.copy()
                for dot in dash.conn_dots:
                    min_dist_ep = find_nearest_geom(dot.point, endpoints)
                    endpoints_wo_dots.remove(min_dist_ep)

                # make virtual dots at the extrapolated location
                # find dashes around the virtual dots
                for ep in endpoints_wo_dots:
                    target_p = get_extrapolated_point(dash_body_line, ep)
                    # plot_line(dash_body_line)
                    # plt.plot(target_p.x, target_p.y, marker="+")
                    dot = Dot(target_p)
                    dash_poly_pairs = dot.search_dash_polygons(polygons_wo_dots_tree, poly_line_dict)
                    if len(dash_poly_pairs) == 0:
                        del Dot.all_dots[id(target_p)]
                    else:
                        dash_pairs = create_dash_pairs(dot, dash_poly_pairs)  # Dash object is created
                        dot.save_dash_pairs(dash_pairs)

    logging.info('Extracting the lines from dots and dash polygons')
    # obtain dash lines
    for dash in Dash.all_dashes.values():
        if id(dash.poly) in poly_line_dict:
            dash_body_line = poly_line_dict[id(dash.poly)]
        else:
            dash_body_line = create_centerline(dash.poly)
            poly_line_dict[id(dash.poly)] = dash_body_line

        if dash_body_line is not None:
            endpoints = get_line_endpoints(dash_body_line)

            # prepare dash body line to be merged
            if isinstance(dash_body_line, MultiLineString):
                lines = list(dash_body_line.geoms)
            else:
                lines = [dash_body_line]

            # make the connecting line from dash's endpoint to dot
            for dot in dash.conn_dots:
                min_dist_ep = find_nearest_geom(dot.point, endpoints)
                lines.append(LineString([min_dist_ep, dot.point]))

            # merged line with the dash body line and connecting line to dots
            mline = linemerge(lines)

            # print(mline)
            # plot_line(mline)
            # plt.show()
            # find shortest paths between dots
            logging.info('Finding shortest path between dots to connect dots and dashes')
            path_lines = [get_path_line(mline, dot1.point, dot2.point) for dot1, dot2 in
                          combinations(dash.conn_dots, 2)]

            if len(path_lines) > 0:
                # TODO: temporary exception handling for debugging purpose
                try:
                    dash_line = unary_union(path_lines)
                    dash.dash_line = dash_line

                    # plot the final line for dash
                    #plot_line(dash_line)

                    if isinstance(dash_line, MultiLineString):
                        all_drawn_lines.extend(list(dash_line.geoms))
                    elif isinstance(dash_line, LineString):  # LineString
                        all_drawn_lines.append(dash_line)
                    else:
                        TypeError(
                            f'Inappropriate type: {type(dash_line)} for dash_line whereas a MultiLineString or LineString is expected')

                except Exception as err:
                    logging.debug(err)
                    for line in path_lines:
                        logging.debug(line)

    logging.info('Merging all extracted lines')
    return linemerge(all_drawn_lines)
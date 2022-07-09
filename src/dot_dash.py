from shapely import affinity
from shapely.geometry import Point, box, Polygon
from settings import DASH_SEARCH_BOX_W, DASH_SEARCH_BOX_H, MAX_DOT_LENGTH
from line_proc import is_endpoint_inside, create_centerline

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
                    if id(poly) in poly_line_dict:
                        poly_centerline = poly_line_dict[id(poly)]
                    else:
                        poly_centerline = create_centerline(poly)
                        poly_line_dict[id(poly)] = poly_centerline

                    if is_endpoint_inside(poly_centerline, sbox) and poly.length > MAX_DOT_LENGTH:
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


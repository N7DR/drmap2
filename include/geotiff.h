// $Id: geotiff.h 2 2019-10-03 18:02:04Z n7dr $

// Source code copyright 2003, 2005 IPfonix, Inc.
// Unauthorized copying strictly prohibited

#ifndef GEOTIFFH
#define GEOTIFFH

#include "macros.h"
#include "string_functions.h"
#include "x_error.h"

#include <array>
#include <atomic>
#include <string>
#include <string_view>

#include <gdal/gdal.h>
#include <gdal/gdal_priv.h>

// error numbers
constexpr int GEOTIFF_NODATA { -1 };

constexpr double RE   { 6371000.0 };                  // radius in m
constexpr double PI   { 3.14159265358979 };
constexpr double DTOR { PI / 180.0 };
constexpr double RTOD { 1.0 / DTOR };

class geotiff_tile;

// forward declarations
void        download_tiff_file_if_necessary(const int llc, const std::string_view local_directory);
std::string local_tile_filename(const int llc, const std::string_view local_directory);

// -----------  geotiff_transform ----------------

/*! \class  geotiff_transform
    \brief  Encapsulate a GeoTransform
*/

class geotiff_transform
{
protected:

/*
  https://gdal.org/en/stable/tutorials/geotransforms_tut.html#
  https://gdal.org/en/stable/tutorials/raster_api_tut.html#opening-the-file

[GT(0) x-coordinate of the upper-left corner of the upper-left pixel.
GT(1) w-e pixel resolution / pixel width.
GT(2) row rotation (typically zero).
GT(3) y-coordinate of the upper-left corner of the upper-left pixel.
GT(4) column rotation (typically zero).
GT(5) n-s pixel resolution / pixel height (negative value for a north-up image).
]
*/
  std::array<double, 6>  _data { };       ///< the actual transformation data

  bool _north_up_image { false };         ///< is this a north-up image?

public:

/// default constructor
  inline geotiff_transform(void)
    { };

/// create from a GDAL GDALDatasetUniquePtr
  explicit geotiff_transform(const GDALDatasetUniquePtr ptr);

  READ(north_up_image);         ///< is this a north-up image?

/// return a pointer to the start of the array
  inline double* data_ptr(void)
    { return static_cast<double*>(_data.data()); }

/// return the x coordinate of the upper-left corner of the upper-left pixel
  inline double x_ul(void) const
    { return _data[0]; }

/// return the width of a pixel
  inline double px_w(void) const
    { return _data[1]; }

/// how much should rows be rotated?
  inline double row_rot(void) const
    { return _data[2]; }

/// return the y coordinate of the upper-left corner of the upper-left pixel
  inline double y_ul(void) const
    { return _data[3]; }

/// how much should columns be rotated?
  inline double column_rot(void) const
    { return _data[4]; }

/// return the height of a pixel (always +ve)
  inline double px_h(void) const
    { return _data[5]; }

/*! \brief  Correct for negative value of pixel height
*/
  void normalise(void);

/// return the object as a human-readable string
  std::string to_string(void) const;
};

// -----------  geotiff_tile ----------------

/*! \class  geotiff_tile
    \brief  Encapsulate a USGS GeoTIFF tile
*/

class geotiff_tile
{ enum class QUADRANT { Q0,     ///< too close to the centre for the quadrant number to be meaningful
                        Q1,     ///< first quadrant
                        Q2,     ///< second quadrant
                        Q3,     ///< third quadrant
                        Q4      ///< fourth quadrant
                      };        ///< quadrants in a tile

protected:

  std::string _data_filename { };           ///< name of the data file
  int         _n_bands       { 0 };         ///< number of bands in the tile
  int         _n_columns     { 0 };         ///< number of columns in the tile
  int         _n_rows        { 0 };         ///< number of rows in the tile
  bool        _valid         { false };     ///< does the tile contain data?

  double _xl { 0 };             ///< longitude of western edge
  double _xr { 0 };             ///< longitude of eastern edge
  double _yb { 0 };             ///< latitude of southern edge
  double _yt { 0 };             ///< latitude of northern edge

// assume float is 32-bit (this is checked before use)
  std::vector<std::vector<float>> _data { };    ///< actual data in the tile; enables access as [latitude][longitude]

  double _no_data { 0 };        ///< the nodata value

  GDALDatasetUniquePtr _data_ptr { nullptr };     ///< pointer to the data in the tile

  geotiff_transform _transform { };               ///< the transform values for the tile

  static std::mutex        _geotiff_mutex;          // used to ensure there is no race condition when checking whether GDAL is registered
  static std::atomic<bool> _registered;             // used to keep track of whether GDAL is registered

/*! \brief          Is a value between two other values?
    \param  value   value to test
    \param  v1      one bound
    \param  v2      other bound
    \return         whether <i>value</i> is between <i>v1</i> and <i>v2</i>

    <i>v1</i> and <i>v2</i> may be in either order.
    Perhaps take this outside the class?
*/
 template<typename T>
   bool _is_between(const T& value, const T& v1, const T& v2) const
   { if ( (value >= v1) and (value <= v2) )
       return true;

     if ( (value >= v2) and (value <= v1) )
       return true;

     return false;
   }

/*! \brief              Map a latitude to a row number
    \param  latitude    latitude to map
    \return             row number that contains latitude <i>latitude</i>
*/
  inline int _map_latitude_to_index(const double latitude) const
    { return static_cast<int>( (_yt - latitude) / px_h() ); }

/*! \brief              Map a longitude to a column number
    \param  longitude   longitude to map
    \return             column number that contains longitude <i>longitude</i>
*/
  inline int _map_longitude_to_index(const double longitude) const
    { return static_cast<int>( (longitude - _xl) / px_w() ); }

/*! \brief              Return the QUADRANT within a cell for a particular latitude and longitude
    \param  latitude    latitude of point
    \param  longitude   longitude of point
    \return             the quadrant within a tile cell in which the point lies

    Returns Q0 if the point is within one metre of the centre of the cell
*/
  QUADRANT _quadrant(const double latitude, const double longitude) const;

/// Return a quadrant as a readable string (e.g., "Q3")
  std::string _quadrant_name(const QUADRANT q) const;

public:

/// default constructor
  inline geotiff_tile(void)
    { }

/// construct from file
  explicit geotiff_tile(const std::string_view fn);

  READ(data_filename)      ///< name of the data file
  READ(n_bands);           ///< number of bands in the tile
  READ(n_columns);         ///< number of columns in the tile
  READ(n_rows);            ///< number of rows in the tile
  READ(no_data);           ///< the nodata value
  READ(valid);             ///< does the tile contain data?

/// Return the filename associated with the tile
  inline std::string filename(void) const
    { return data_filename(); }

/// Return the name of the driver associated with the tile
  inline std::string driver_name(void) const
    { return ( _valid? std::string(_data_ptr -> GetDriverName()) : std::string { } ); }     // https://gdal.org/en/stable/doxygen/classGDALDataset.html

/// return the x coordinate of the upper-left corner of the upper-left pixel
  inline double x_ul(void) const
    { return _transform.x_ul(); }

/// return the width of a pixel
  inline double px_w(void) const
    { return _transform.px_w(); }

/// how much should rows be rotated?
  inline double row_rot(void) const
    { return _transform.row_rot(); }

/// return the y coordinate of the upper-left corner of the upper-left pixel
  inline double y_ul(void) const
    { return _transform.y_ul(); }

/// how much should columns be rotated?
  inline double column_rot(void) const
    { return _transform.column_rot(); }

/// return the height of a pixel (always +ve)
  inline double px_h(void) const
    { return _transform.px_h(); }

/*! \brief             Convert a latitude and longitude to the equivalent indices
     \param  latitude   latitude of point
     \param  longitude  longitude of point
     \return            indices equivalent to <i>latitude</i> and <i>longitude</i>
*/
  inline std::pair<int, int> index_pair(const double latitude, const double longitude) const
    { return { _map_latitude_to_index(latitude), _map_longitude_to_index(longitude) }; }

/*! \brief          Get the central latitude and longitude for a cell identified by an index pair
    \param  ipair   index pair (latitude index, longitude index) -- indices are wrt 0
    \return         the latitude and longitude of the centre of the cell [ipair.first][ipair.second]

    Performs no bounds checking
*/
  std::pair<double /* latitude */, double /* longitude */> cell_centre(const std::pair<int, int>& ipair) const;

/*! \brief          Get the central latitude and longitude for a cell identified by an index pair
    \param  ipair   index pair (latitude index, longitude index)
    \return         the latitude and longitude of the centre of the cell [ipair.first][ipair.second]
*/
  inline std::pair<double, double> cell_centre(const std::pair<double, double>& ipair) const
    { return cell_centre(index_pair(ipair.first, ipair.second)); }

/*! \brief      Get the central latitude and longitude for a cell containing a particular latitude and longitude
    \param      latitude   latitude of point
    \param      longitude  longitude of point    \param  ipair   index pair (latitude index, longitude index)
    \return     the latitude and longitude of the centre of the cell that contains [latitude, longitude]
*/
  inline std::pair<double, double> cell_centre(const double latitude, const double longitude) const
    { return cell_centre(index_pair(latitude, longitude)); }

/*! \brief              The value of the cell that contains a point
    \param  latitude    latitude of point
    \param  longitude   longitude of point
    \return             the value of the cell containing the point

    Returns <i>_no_data</i> if the point is not within the tile
*/
  float cell_value(const double latitude, const double longitude) const;

/*! \brief      The value of the cell with particular indices
    \param  ip  index pair (latitude index, longitude index)
    \return     the value of the cell [ip.first][ip.second]

    Performs no bounds checking
*/
  inline float cell_value(const std::pair<int, int>& ip) const  // pair is lat index, long index
    { return _data[ip.first][ip.second]; }

/*! \brief              The weighted average of the cells near a point
    \param  latitude    latitude of point
    \param  longitude   longitude of point
    \return             the interpolated value of nearby cells, to give the value for the point

    Returns <i>_no_data</i> if the point is not within the tile, or if there are insufficient data
    to return a valid response.
*/
  float interpolated_value(const double latitude, const double longitude) const;

/*! \brief      The weighted average of the cells near a point
    \param  ll  latitude and longitude of point
    \return     the interpolated value of nearby cells, to give the value for the point

    Returns <i>_no_data</i> if the point is not within the tile
*/
  inline float interpolated_value(const std::pair<double, double>& ll) const
    { return interpolated_value(ll.first, ll.second); }

/*! \brief              Is a point within the tile?
    \param  latitude    latitude of point
    \param  longitude   longitude of point
    \return             whether the point is within the tile
*/
  inline bool is_in_tile(const double latitude, const double longitude) const
    { return (_is_between(latitude, _yb, _yt) and _is_between(longitude, _xl, _xr)); }

/// does a height appear to be valid?
  inline bool valid_height(const float h) const
    { return ( h > (-9999 + 1) ); }

  friend class geotiff_band;
};

// -----------  geotiff_band ----------------

/*! \class  geotiff_band
    \brief  Encapsulate a GeoTIFF raster band
*/

class geotiff_band
{
protected:

  GDALRasterBand* _data_ptr { nullptr };    // pointer to the data in the band

  double _no_data { 0 };        ///< the default nodata value


public:

/*! \brief            Constructor
    \param  gt        a tile
    \param  band_nr   the number of the band in the tile, wrt 1
*/
  inline geotiff_band(const geotiff_tile& gt, const int band_nr /* wrt 1 */) :
    _data_ptr(gt._data_ptr -> GetRasterBand(band_nr)),
    _no_data(_data_ptr -> GetNoDataValue())
  { }
};

/*! \brief              Obtain the bearing (from north) associated with displacement by an amount horizontally and vertically
    \param  delta_x     number and direction of horizontal units
    \param  delta_y     number and direction of vertical units
    \return             the bearing, in degrees, associated with displacement by <i>delta_x</i> and <i>delta_y</i>

    "bearing" here means the initial bearing at which one must leave the central point; it is also the bearing at which a cell
    corresponding to <i>delta_x</i> and <i>delta_y</i> is to be plotted
*/
double bearing(const int delta_x, const int delta_y);  // bearing in degrees

/*! \brief              Obtain latitude and longitude corresponding to a bearing and distance from a point
    \param  lat1        latitude of source, in degrees (+ve north)
    \param  long1       longitude of source, in degrees (+ve east)
    \param  bearing_d   bearing in degrees from source
    \param  distance_m  distance in metres (along Earth's surface) from source
    \return             latitude and longitude of target

Formula: φ2 = asin( sin φ1 ⋅ cos δ + cos φ1 ⋅ sin δ ⋅ cos θ )
	λ2 = λ1 + atan2( sin θ ⋅ sin δ ⋅ cos φ1, cos δ − sin φ1 ⋅ sin φ2 )
where 	φ is latitude, λ is longitude, θ is the bearing (clockwise from north), δ is the angular distance d/R; d being the distance travelled, R the earth’s radius
JavaScript:
(all angles
in radians)

var φ2 = Math.asin( Math.sin(φ1)*Math.cos(d/R) +
                    Math.cos(φ1)*Math.sin(d/R)*Math.cos(brng) );
var λ2 = λ1 + Math.atan2(Math.sin(brng)*Math.sin(d/R)*Math.cos(φ1),
                         Math.cos(d/R)-Math.sin(φ1)*Math.sin(φ2));

	The longitude can be normalised to −180…+180 using (lon+540)%360-180
Excel:
(all angles
in radians)
	lat2: =ASIN(SIN(lat1)*COS(d/R) + COS(lat1)*SIN(d/R)*COS(brng))
lon2: =lon1 + ATAN2(COS(d/R)-SIN(lat1)*SIN(lat2), SIN(brng)*SIN(d/R)*COS(lat1))
*/
std::pair<double, double> ll_from_bd(const double lat1, const double long1, const double bearing_d, const double distance_m);

/*! \brief              Obtain latitude and longitude corresponding to a bearing and distance from a point
    \param  ll          latitude and longitude of source, in degrees (+ve north, +ve east)
    \param  bearing_d   bearing in degrees from source
    \param  distance_m  distance in metres (along Earth's surface) from source
    \return             latitude and longitude of target
*/
inline std::pair<double, double> ll_from_bd(const std::pair<double, double>& ll, const double& bearing_d /* degrees */, const double& distance_m /* metres */)
  { return ll_from_bd(ll.first, ll.second, bearing_d, distance_m); }

/*! \brief          Obtain distance in km between two locations
    \param  lat1    latitude of source, in degrees (+ve north)
    \param  long1   longitude of source, in degrees (+ve east)
    \param  lat2    latitude of target, in degrees (+ve north)
    \param  long2   longitude of target, in degrees (+ve east)
    \return         distance between source and target, in km

    See http://www.movable-type.co.uk/scripts/latlong.html:

    a = sin²(Δφ/2) + cos(φ1).cos(φ2).sin²(Δλ/2)
    c = 2.atan2(√a, √(1−a))
    d = R.c
    where   φ is latitude, λ is longitude, R is earth’s radius (mean radius = 6,371km)

    θ = atan2( sin(Δλ).cos(φ2), cos(φ1).sin(φ2) − sin(φ1).cos(φ2).cos(Δλ) )
*/
double distance(const double lat1, const double long1, const double lat2, const double long2);

/*! \brief                Obtain distance in km between two locations
    \param  lat_long_1    latitude and longitude of first point, in degrees (+ve north, +ve east)
    \param  lat_long_2    latitude and longitude of second point, in degrees (+ve north, +ve east)
    \return               distance between the first and second point, in km
*/
inline double distance(const std::pair<double, double>& lat_long_1, const std::pair<double, double>& lat_long_2)
  { return distance(lat_long_1.first, lat_long_1.second, lat_long_2.first, lat_long_2.second); }


/*! \brief              Obtain distance in km between two locations
    \param  lat_long_1  latitude and longitude of first point, in degrees (+ve north, +ve east)
    \param  lat2        latitude of second point, in degrees (+ve north)
    \param  long2       longitude of second point, in degrees (+ve east)
    \return             distance between the first and second point, in km
*/
inline double distance(const std::pair<double, double>& lat_long_1, const double lat2, const double long2)
  { return distance(lat_long_1.first, lat_long_1.second, lat2, long2); }

/*! \brief              Obtain distance in km between two locations
    \param  lat1        latitude of first point, in degrees (+ve north)
    \param  long1       longitude of first point, in degrees (+ve east)
    \param  lat_long_2  latitude and longitude of second point, in degrees (+ve north, +ve east)
    \return             distance between the first and second point, in km
*/
inline double distance(const double lat1, const double long1, const std::pair<double, double>& lat_long_2)
  { return distance(lat1, long1, lat_long_2.first, lat_long_2.second); }

/*! \brief              Return a base filename derived from latitude and longitude
    \param  latitude    latitude
    \param  longitude   longitude
    \return             the base name of the file that contains the data for the tile that contains the point at <i>latitude</i>, <i>longitude</i>

    "nLLwLLL"
*/
std::string base_filename(const double latitude, const double longitude);

// "nLLwLLL"
inline std::string base_filename(const std::pair<double, double>& ll)
  { return base_filename(ll.first, ll.second); }

/*! \brief              Return a base filename derived from latitude and longitude
    \param  llcode      the llcode [lat * 1000 + (+ve)long]
    \return             the base name of the file that contains the data for the tile that contains the point at <i>latitude</i>, <i>longitude</i>

    "nLLwLLL"
*/
std::string base_filename(const int llcode);

// lambdas can't be overloaded (as they have no user-accessible name), so we have to define functions to convert various things to lat-long code

/*! \brief              Obtain a lat-long code from latitude and longitude
    \param  latitude    latitude
    \param  longitude   longitude
    \return             llcode [lat * 1000 + (+ve)long]
*/
inline int llc(const double latitude, const double longitude)
  { return ( int(latitude + 1) * 1000 + int(-(longitude - 1) ) ); }

/*! \brief      Obtain a lat-long code from latitude and longitude
    \param  ll  latitude and longitude
    \return     llcode [lat * 1000 + (+ve)long]
*/
inline int llc(const std::pair<double, double>& ll)
  { return ( llc(ll.first, ll.second) ); };

/*! \brief                Obtain a lat-long code from a filename
    \param  basefilename  filename in the form "nLLwLLL"
    \return               llcode [lat * 1000 + (+ve)long]
*/
inline int llc(const std::string_view basefilename) //"nLLwLLL"
  { return ( from_string<int>(basefilename.substr(1, 2)) * 1000 + from_string<int>(basefilename.substr(4, 3))); }

/*! \brief      Correction to be applied to a distace due to the curvature of the Earth
    \param  d   original distance
    \return     amount by which <i>d</i> has to be corrected
*/
inline double curvature_correction(const double d)
  { return ( ( 1 - cos(d / RE) ) * RE ); }

// -----------  geotiff_error ----------------

/*! \class  geotiff_error
    \brief  Errors related to GeoTIFF objects
*/
class geotiff_error : public x_error
{
protected:

public:

/*!	\brief	    Construct from error code and reason
	\param	n	error code
	\param	s	reason
*/
  geotiff_error(const int n, const std::string& s) :
    x_error(n, s)
  { }
};

#endif    // GEOTIFFH


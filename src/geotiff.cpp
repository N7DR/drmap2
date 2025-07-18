// $Id: geotiff.cpp 2 2019-10-03 18:02:04Z n7dr $

// Source code copyright 2003, 2005 IPfonix, Inc.
// Unauthorized copying strictly prohibited

#include "diskfile.h"
#include "geotiff.h"
#include "string_functions.h"

#include <iostream>
#include <thread>

#include <gdal/gdal.h>

using namespace std;
using namespace std::this_thread;

extern bool debug;

// globals
mutex        geotiff_tile::_geotiff_mutex { };          // used to ensure there is no race condition when checking whether GDAL is registered
atomic<bool> geotiff_tile::_registered    { false };    // used to keep track of whether GDAL is registered

/*! \brief                      Get the filename for a local GeoTIFF tile
    \param  llc                 the llcode [lat * 1000 + (+ve)long]
    \param  local_directory     the local directory containing USGS files
*/
string local_tile_filename(const int llc, const string_view local_directory)
{ const string nwcode            { base_filename(llc) };                      // something like: n40w105
  const string base_name         { "USGS_13_" + nwcode + ".tif"s };           // base name of the geotiff file, including "extension"
  const string local_dirname     { dirname_with_slash(local_directory) };     // ensure a terminating slash
  const string full_geotiff_name { local_dirname + base_name };               // full name of local geotiff file

  return full_geotiff_name;
}

/*! \brief                      Download a tile from the USGS if we don't already have it
    \param  llc                 the llcode [lat * 1000 + (+ve)long]
    \param  local_directory     the local directory containing USGS files
*/
void download_tiff_file_if_necessary(const int llc, const string_view local_directory)
{ constexpr std::chrono::duration DOWNLOAD_RETRY_TIME { 5s };         // time to wait between attempts to download

  bool need_to_download { false };

  const string nwcode            { base_filename(llc) };                              // something like: n40w105
  const string base_name         { "USGS_13_" + nwcode + ".tif"s };                   // base name of the geotiff file, including "extension"
  const string local_dirname     { dirname_with_slash(local_directory) };     // ensure a terminating slash
  const string full_geotiff_name { local_dirname + base_name };               // full name of local geotiff file

  if (!file_exists(full_geotiff_name) or file_empty(full_geotiff_name))
    need_to_download = true;                                                                                    // don't need to download if the header and data files are present

  if (!need_to_download)
    return;                                                                                                     // we're done here

  directory_create_if_necessary(local_directory);

  const string remote_directory { R"(https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/13/TIFF/current/)" + nwcode + "/"s };
  const string remote_filename  { base_name };
  const string local_filename   { full_geotiff_name };

  bool   downloaded { false };
  string command    { };

  constexpr int MAX_TRIES { 10 };   // maximum number of times to attempt to download

  for (unsigned int n { 0 }; !downloaded and (n < MAX_TRIES); ++n)
  { if (!file_exists(local_filename))                                          // don't download if it already exists
    { cout << "File " << local_filename << " does not exist; download attempt number " << (n + 1) << endl;

      command = R"(wget -q -O )"  + local_filename + " "s + remote_directory + remote_filename;

      if (debug)
        cout << "command = ***" << command << "***" << endl;

      system(command.c_str());

      if (debug)
        cout << "command completed: " << command.c_str() << endl;
    }
    else
    { if (debug)
        cout << "file already exists: " << local_filename << endl;
    }

// did we succeed?
    if (file_exists(local_filename))
    { if (file_empty(local_filename))
        file_delete(local_filename);
      else                              // exists and not empty; assume download was OK
        downloaded = true;
    }

    if (!downloaded)
      sleep_for(DOWNLOAD_RETRY_TIME);          // USGS seems to need requests slowed down, otherwise they all instantly fail
  }

  cout << (downloaded ? "Download succeeded" : "Download did not succeed") << endl;

  if (!downloaded)
  { cerr << "download failed; exiting" << endl;
    exit(-1);
  }
}

// -----------  geotiff_transform ----------------

/*! \class  geotiff_transform
    \brief  Encapsulate a GeoTransform
*/

/// create from a GDAL GDALDatasetUniquePtr
geotiff_transform::geotiff_transform(const GDALDatasetUniquePtr ptr)
{ const int status { ptr -> GetGeoTransform(data_ptr()) };

  if (status == CE_Failure)
  { cerr << "unable to obtain a GeoTransform array" << endl;
    exit(-1);
  }

  normalise();
}

/// return the object as a human-readable string
string geotiff_transform::to_string(void) const
{ string rv;

  rv += "[0]:  x-coordinate of the upper-left corner of the upper-left pixel = " + ::to_string(_data[0]) + EOL;
  rv += "[1]:  w-e pixel resolution / pixel width                            = " + ::to_string(_data[1]) + EOL;
  rv += "[2]:  row rotation                                                  = " + ::to_string(_data[2]) + EOL;
  rv += "[3]:  y-coordinate of the upper-left corner of the upper-left pixel = " + ::to_string(_data[3]) + EOL;
  rv += "[4]:  column rotation                                               = " + ::to_string(_data[4]) + EOL;
  rv += "[5]:  n-s pixel resolution / pixel height                           = " + ::to_string(_data[5]) + EOL;
  rv += "north-up image                                                      =  " + (_north_up_image ? "TRUE"s : "FALSE"s);

  return rv;
}

/*! \brief  Correct for negative value of pixel height
*/
void geotiff_transform::normalise(void)
{ if (_data[5] < 0)
  { _data[5] = -(_data[5]);
    _north_up_image = true;
  }
}

// -----------  geotiff_tile ----------------

/*! \class  geotiff_tile
    \brief  Encapsulate a USGS GeoTIFF tile
*/

/*! \brief              Return the QUADRANT within a cell for a particular latitude and longitude
    \param  latitude    latitude of point
    \param  longitude   longitude of point
    \return             the quadrant within a tile cell in which the point lies

    Returns Q0 if the point is within one metre of the centre of the cell

    NB: longitudes are generally -ve; more -ve => more west

    Q1: NE
    Q2: NW
    Q3: SW
    Q4: SE
*/
geotiff_tile::QUADRANT geotiff_tile::_quadrant(const double latitude, const double longitude) const
{ const auto   centre                    { cell_centre(latitude, longitude) };
  const auto&  [c_latitude, c_longitude] { centre };
  const double dist                      { distance(centre, latitude, longitude) };

// arbitrary, if distance < 1m, then treat as 0
  if (dist < 1)
    return QUADRANT::Q0;

  if ( (latitude >= c_latitude) and (longitude >= c_longitude) )
    return QUADRANT::Q1;

  if ( (latitude >= c_latitude) and (longitude <= c_longitude) )
    return QUADRANT::Q2;

  if ( (latitude <= c_latitude) and (longitude <= c_longitude) )
    return QUADRANT::Q3;

  return QUADRANT::Q4;
}

/// Return a quadrant as a readable string (e.g., "Q3")
string geotiff_tile::_quadrant_name(const geotiff_tile::QUADRANT q) const
{ switch (q)
  { case QUADRANT::Q0 : return "Q0: near centre"s;
    case QUADRANT::Q1 : return "Q1: NE"s;
    case QUADRANT::Q2 : return "Q2: NW"s;
    case QUADRANT::Q3 : return "Q3: SW"s;
    case QUADRANT::Q4 : return "Q4: SE"s;
  }

// keep stupid compiler happy
  return "ERROR"s;
}

/// construct from file
geotiff_tile::geotiff_tile(const string_view fn) :
  _data_filename(fn)
{
// check that GDAL is ready to use
  { lock_guard lg(_geotiff_mutex);

    if (!_registered)
    { GDALAllRegister();
      _registered = true;
    }
  }

  _data_ptr = GDALDatasetUniquePtr(GDALDataset::FromHandle(GDALOpen( string(fn).c_str(), GA_ReadOnly )));

  if (_data_ptr)
    _valid = true;
  else
  { cerr << "Unable to obtain GDAL data from file: " << fn << endl;
    exit(-1);
  }

  _n_bands = _data_ptr -> GetRasterCount();

  if (_n_bands != 1)                                                        // check that we have only one band
  { cerr << "Invalid number of bands in tile: " << _n_bands << endl;
    exit(-1);
  }

  _n_columns = _data_ptr -> GetRasterXSize();
  _n_rows    = _data_ptr -> GetRasterYSize();

  int status { _data_ptr -> GetGeoTransform(_transform.data_ptr()) };   // fill _transform data (stupid C-style function that takes a destination address)

  if (status)
  { _valid = false;
    return;
  }

  _transform.normalise();

  _xl = x_ul();
  _xr = _xl + px_w() * n_columns();
  _yt = y_ul();
  _yb = _yt - px_h() * n_rows();

// create the space needed to hold the data
  _data.reserve(_n_rows);

  for (int n { 0 }; n < _n_rows; ++n)
  { _data += vector<float> { };

    _data[n].reserve(_n_columns);
  }

  float* f_ptr = new float [_n_columns];

  GDALRasterBand* poBand { _data_ptr -> GetRasterBand( 1 ) };

  _no_data = poBand -> GetNoDataValue(&status);

  if (!status)    // failure
  { _valid = false;
    _no_data = 0;
  }

/*
 CPLErr GDALRasterBand::RasterIO( GDALRWFlag eRWFlag,
                                int nXOff, int nYOff, int nXSize, int nYSize,
                                void * pData, int nBufXSize, int nBufYSize,
                                GDALDataType eBufType,
                                int nPixelSpace,
                                int nLineSpace )
*/

  if (debug)
    cout << "tile file: " << filename() << endl;

  for (int n { 0 }; n < _n_rows; ++n)
  { const int status { poBand -> RasterIO( GF_Read, 0, n, _n_columns, 1, f_ptr, _n_columns, 1, GDT_Float32, 0, 0 ) };

    if (status)
    { cout << "error reading line from tile file; status = " << status << endl;
      exit(-1);
    }

    _data[n].assign(f_ptr, f_ptr + _n_columns);  // ugh!
  }

  if (debug)
    cout << "data filled for tile filename: " << fn << "; data[0][0] = " << _data[0][0] << endl;

  delete [] f_ptr;
}

/*! \brief          The latitude and longitude of the cell with particular indices
    \param  ipair   index pair (latitude index, longitude index) -- indices are wrt 0
    \return         the latitude and longitude of the centre of the cell [ipair.first][ipair.second]

    Performs no bounds checking
*/
pair<double, double> geotiff_tile::cell_centre(const pair<int, int>& ipair) const
{ const auto& [lat_index, lon_index] { ipair };

  const double lat_0     { _yt - px_h() / 2 };
  const double latitude  { lat_0 - (lat_index * px_h()) };
//  const double latitude  { lat_0 - ipair.first * px_h() };
  const double long_0    { _xl + px_w() / 2 };
 // const double longitude { long_0 + ipair.second * px_w() };
  const double longitude { long_0 + (lon_index * px_w()) };

  return { latitude, longitude };
}

/*! \brief              The value of the cell that contains a point
    \param  latitude    latitude of point
    \param  longitude   longitude of point
    \return             the value of the cell containing the point

    Returns <i>_no_data</i> if the point is not within the tile
*/
float geotiff_tile::cell_value(const double latitude, const double longitude) const
{ if (is_in_tile(latitude, longitude))
  { const int row_nr    { _map_latitude_to_index(latitude) };
    const int column_nr { _map_longitude_to_index(longitude) };

    return _data[row_nr][column_nr];
  }

  return no_data();
}

/*! \brief      The value of the cell with particular indices
    \param  ip  index pair (latitude index, longitude index)
    \return     the value of the cell [ip.first][ip.second]

    Performs no bounds checking
*/
//float geotiff_tile::cell_value(const std::pair<int, int>& ip) const  // pair is lat index, long index
//{ return _data[ip.first][ip.second];
//}

/*! \brief              The weighted average of the cells near a point
    \param  latitude    latitude of point
    \param  latitude    latitude of point
    \param  longitude   longitude of point
    \return             the interpolated value of nearby cells, to give the value for the point

    Returns <i>_no_data</i> if the point is not within the tile, or if there are insufficient data
    to return a valid response.
*/
float geotiff_tile::interpolated_value(const double latitude, const double longitude) const
{ const pair<double, double> ll_centre { cell_centre(latitude, longitude) };
  const QUADRANT             q         { _quadrant(latitude, longitude) };

  switch (q)
  { case QUADRANT::Q0 :
    { const auto cv { cell_value(latitude, longitude) };

      if (cv < -9000)
        throw geotiff_error(GEOTIFF_NODATA, ( "Q0: Insufficient data when interpolating at "s + ::to_string(latitude) + ", "s + ::to_string(longitude)) );

      return cv;
    }

    case QUADRANT::Q1 :
    { array<float, 4>  corner_heights;
      array<double, 4> corner_distances;

      const pair<int, int> indices_0 { index_pair(latitude, longitude) };

      corner_heights[0] = cell_value(latitude, longitude);             // height at centre
      corner_distances[0] = distance(ll_centre, latitude, longitude);  // distance from centre to lat/lon

// move right => increment second index
      const pair<int, int>       indices_1   { indices_0.first, indices_0.second + 1 };
      const pair<double, double> ll_centre_1 { cell_centre(indices_1) };

      corner_heights[1] = cell_value(indices_1);
      corner_distances[1] = distance(ll_centre_1, latitude, longitude);

// move up => decrement first index
      const pair<int, int>       indices_2   { indices_0.first - 1, indices_0.second };
      const pair<double, double> ll_centre_2 { cell_centre(indices_2) };

      corner_heights[2] = cell_value(indices_2);
      corner_distances[2] = distance(ll_centre_2, latitude, longitude);

// move up and to the right => decrement first index, increment second index
      const pair<int, int>       indices_3   { indices_0.first - 1, indices_0.second + 1 };
      const pair<double, double> ll_centre_3 { cell_centre(indices_3) };

      corner_heights[3] = cell_value(indices_3);
      corner_distances[3] = distance(ll_centre_3, latitude, longitude);

// calculate the weighted mean
      double h     { 0 };
      double inv_d { 0 };

      for (auto n = 0; n < 4; ++n)
      { if (valid_height(corner_heights[n]))
        { h += corner_heights[n] / corner_distances[n];
          inv_d += 1 / corner_distances[n];
        }
      }

      if (inv_d == 0)
        throw geotiff_error(GEOTIFF_NODATA, ( "Q1: Insufficient data when interpolating at "s + ::to_string(latitude) + ", " + ::to_string(longitude)) );

      const double rv { (inv_d == 0) ? -9999 : (h / inv_d) };

      return rv;
    }

    case QUADRANT::Q2 :
    { array<float, 4> corner_heights;
      array<double, 4> corner_distances;

      const pair<int, int> indices_0 { index_pair(latitude, longitude) };

      corner_heights[0] = cell_value(latitude, longitude);             // height at centre
      corner_distances[0] = distance(ll_centre, latitude, longitude);  // distance from centre to lat/lon

// move left => decrement second index
      const pair<int, int>       indices_1   { indices_0.first, indices_0.second - 1 };
      const pair<double, double> ll_centre_1 { cell_centre(indices_1) };

      corner_heights[1] = cell_value(indices_1);
      corner_distances[1] = distance(ll_centre_1, latitude, longitude);

// move up => decrement first index
      const pair<int, int>       indices_2   { indices_0.first - 1, indices_0.second };
      const pair<double, double> ll_centre_2 { cell_centre(indices_2) };

      corner_heights[2] = cell_value(indices_2);
      corner_distances[2] = distance(ll_centre_2, latitude, longitude);

// move up and to the left => decrement first index, decrement second index
      const pair<int, int>       indices_3   { indices_0.first - 1, indices_0.second - 1 };
      const pair<double, double> ll_centre_3 { cell_centre(indices_3) };

      corner_heights[3] = cell_value(indices_3);
      corner_distances[3] = distance(ll_centre_3, latitude, longitude);

// calculate the weighted mean
      double h     { 0 };
      double inv_d { 0 };

      for (auto n = 0; n < 4; ++n)
      { h += corner_heights[n] / corner_distances[n];
        inv_d += 1 / corner_distances[n];
      }

      if (inv_d == 0)
        throw geotiff_error(GEOTIFF_NODATA, ( "Q2: Insufficient data when interpolating at "s + ::to_string(latitude) + ", " + ::to_string(longitude)) );

      const double rv { h / inv_d };

      return rv;
    }

    case QUADRANT::Q3 :
    { array<float, 4>  corner_heights;
      array<double, 4> corner_distances;

      const pair<int, int> indices_0 { index_pair(latitude, longitude) };

      corner_heights[0] = cell_value(latitude, longitude);             // height at centre
      corner_distances[0] = distance(ll_centre, latitude, longitude);  // distance from centre to lat/lon

// move left => decrement second index
      const pair<int, int>       indices_1   { indices_0.first, indices_0.second - 1 };
      const pair<double, double> ll_centre_1 { cell_centre(indices_1) };

      corner_heights[1] = cell_value(indices_1);
      corner_distances[1] = distance(ll_centre_1, latitude, longitude);

// move down => increment first index
      const pair<int, int>       indices_2   { indices_0.first + 1, indices_0.second };
      const pair<double, double> ll_centre_2 { cell_centre(indices_2) };

      corner_heights[2] = cell_value(indices_2);
      corner_distances[2] = distance(ll_centre_2, latitude, longitude);

// move down and to the left => increment first index, decrement second index
      const pair<int, int>       indices_3   { indices_0.first + 1, indices_0.second - 1 };
      const pair<double, double> ll_centre_3 { cell_centre(indices_3) };

      corner_heights[3] = cell_value(indices_3);
      corner_distances[3] = distance(ll_centre_3, latitude, longitude);

// calculate the weighted mean
      double h     { 0 };
      double inv_d { 0 };

      for (auto n = 0; n < 4; ++n)
      { h += corner_heights[n] / corner_distances[n];
        inv_d += 1 / corner_distances[n];
      }

      if (inv_d == 0)
        throw geotiff_error(GEOTIFF_NODATA, ( "Q3: Insufficient data when interpolating at "s + ::to_string(latitude) + ", " + ::to_string(longitude)) );

      const double rv { h / inv_d };

      return rv;
    }

    case QUADRANT::Q4 :
    { array<float, 4>  corner_heights;
      array<double, 4> corner_distances;

      const pair<int, int> indices_0 { index_pair(latitude, longitude) };

      corner_heights[0] = cell_value(latitude, longitude);             // height at centre
      corner_distances[0] = distance(ll_centre, latitude, longitude);  // distance from centre to lat/lon

// move right => increment second index
      const pair<int, int>       indices_1   { indices_0.first, indices_0.second + 1 };
      const pair<double, double> ll_centre_1 { cell_centre(indices_1) };

      corner_heights[1] = cell_value(indices_1);
      corner_distances[1] = distance(ll_centre_1, latitude, longitude);

// move down => increment first index
      const pair<int, int>       indices_2   { indices_0.first + 1, indices_0.second };
      const pair<double, double> ll_centre_2 { cell_centre(indices_2) };

      corner_heights[2] = cell_value(indices_2);
      corner_distances[2] = distance(ll_centre_2, latitude, longitude);

// move down and to the right => increment first index, increment second index
      const pair<int, int>       indices_3   { indices_0.first + 1, indices_0.second + 1 };
      const pair<double, double> ll_centre_3 { cell_centre(indices_3) };

      corner_heights[3] = cell_value(indices_3);
      corner_distances[3] = distance(ll_centre_3, latitude, longitude);

// calculate the weighted mean
      double h     { 0 };
      double inv_d { 0 };

      for (auto n = 0; n < 4; ++n)
      { h += corner_heights[n] / corner_distances[n];
        inv_d += 1 / corner_distances[n];
      }

      if (inv_d == 0)
        throw geotiff_error(GEOTIFF_NODATA, ( "Q4: Insufficient data when interpolating at "s + ::to_string(latitude) + ", " + ::to_string(longitude)) );

      const double rv { h / inv_d };

      return rv;
    }
  }

  return no_data();    // just to keep the compiler happy
}

/*! \brief              Obtain the bearing (from north) associated with displacement by an amount horizontally and vertically
    \param  delta_x     number and direction of horizontal units
    \param  delta_y     number and direction of vertical units
    \return             the bearing, in degrees, associated with displacement by <i>delta_x</i> and <i>delta_y</i>

    "bearing" here means the initial bearing at which one must leave the central point; it is also the bearing at which a cell
    corresponding to <i>delta_x</i> and <i>delta_y</i> is to be plotted
*/
double bearing(const int delta_x, const int delta_y)  // bearing in degrees
{ if ( (delta_x == 0) and (delta_y == 0) )
    return 0;

  if (delta_x == 0)
  { if (delta_y > 0)
      return 0;
    return 180;
  }

  if (delta_y == 0)
  { if (delta_x < 0)
      return 270;
    return 90;
  }

  if ( (delta_x > 0) and (delta_y > 0) )    // 0 -- 90
  { if (delta_y > delta_x)                  // 0 -- 45
      return atan((float)delta_x / (float)delta_y) * RTOD;
    return 90 - ( atan((float)delta_y / (float)delta_x) * RTOD );   // 45 -- 90
  }

  if ( (delta_x > 0) and (delta_y < 0) )    // 90 -- 180
  { if (abs(delta_x) > abs(delta_y))         // 90 -- 135
      return 90 + ( atan((float)abs(delta_y) / (float)abs(delta_x)) * RTOD );
    return 180 - ( atan((float)abs(delta_x) / (float)abs(delta_y)) * RTOD );    // 135 -- 180
  }

  if ( (delta_x < 0) and (delta_y < 0) )    // 180 -- 270
  { if (abs(delta_x) < abs(delta_y))         // 180 -- 225
      return 180 + ( atan((float)abs(delta_x) / (float)abs(delta_y)) * RTOD );
    return 270 - ( atan((float)abs(delta_y) / (float)abs(delta_x)) * RTOD );    // 225 -- 270
  }

// only 270 -- 0 remain
  if (abs(delta_x) > abs(delta_y))         // 270 -- 315
      return 270 + ( atan((float)abs(delta_y) / (float)abs(delta_x)) * RTOD );
  return 360 - ( atan((float)abs(delta_x) / (float)abs(delta_y)) * RTOD );  // 315 -- 360
}

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
pair<double, double> ll_from_bd(const double lat1 /* deg */, const double long1 /* deg */ , const double bearing_d /* degrees clockwise from north */, const double distance_m /* metres */)
{ const double delta   { distance_m / RE };
  const double lat1_r  { lat1 * DTOR };
  const double long1_r { long1 * DTOR };
  const double theta   { bearing_d * DTOR };
  const double lat2_r  { asin ( sin (lat1_r * cos(delta) + cos(lat1_r) * sin(delta) * cos(theta)) ) };
  const double long2_r { long1_r + atan2( sin(theta) * sin(delta) * cos(lat1_r), cos(delta) - sin(lat1_r) * sin(lat2_r)) };
  const double lat2_d  { lat2_r * RTOD };
  const double long2_d { long2_r * RTOD };

  return { lat2_d, long2_d };
}

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
double distance(const double lat1, const double long1, const double lat2, const double long2)
{ const double delta_phi   { lat2 - lat1 };
  const double delta_phi_2 { delta_phi / 2 };

  const double delta_lambda   { long2 - long1 };
  const double delta_lambda_2 { delta_lambda / 2 };

  const double a { sin(delta_phi_2 * DTOR) * sin(delta_phi_2 * DTOR) +
                   cos(lat1 * DTOR) * cos(lat2 * DTOR) * sin(delta_lambda_2 * DTOR) * sin(delta_lambda_2 * DTOR)
                 };

  const double c { 2 * atan2(sqrt(a), sqrt(1 - a)) };
  const double d { RE * c };

  return d;
}

/*! \brief              Return a base filename derived from latitude and longitude
    \param  latitude    latitude
    \param  longitude   longitude
    \return             the base name of the file that contains the data for the tile that contains the point at <i>latitude</i>, <i>longitude</i>

    "nLLwLLL"
*/
string base_filename(const double latitude, const double longitude)
{ const string lat_string  { to_string(static_cast<int>(latitude + 1.0)) };                                  // assumes two digits
  const string long_string { pad_string(to_string(static_cast<int>(-longitude + 1.0)), 3, PAD::LEFT, '0') };

  return ("n"s + lat_string + "w"s + long_string);
}

/*! \brief              Return a base filename derived from latitude and longitude
    \param  llcode      the llcode [lat * 1000 + (+ve)long]
    \return             the base name of the file that contains the data for the tile that contains the point at <i>latitude</i>, <i>longitude</i>

    "nLLwLLL"
*/
std::string base_filename(const int llcode)
{ const string lat_string  { pad_string(to_string(llcode / 1000), 2, PAD::LEFT, '0') };                                  // assumes two digits
  const string long_string { pad_string(to_string(llcode % 1000), 3, PAD::LEFT, '0') };

  return ("n"s + lat_string + "w"s + long_string);
}

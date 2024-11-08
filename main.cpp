
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <limits>
#include <vector>

#include <syncstream>
#define thdebug(x) std::osyncstream(std::cerr) << __FILE__ << " (" << __LINE__ << ") " << #x << " == " << (x) << std::endl;
#define thdebugid(id,x) std::osyncstream(std::cerr) << __FILE__ << " [" << id << "] (" << __LINE__ << ") " << #x << " == " << (x) << std::endl;


#define debug(x)   std::cout << "(" << __LINE__ << ") " << #x << " == " << x << std::endl;

template <typename T>
std::ostream & operator << (std::ostream & sos, const std::vector<T> & data)
{
   sos << "[ ";
   for (const auto & d : data)
      sos << d << " ";
   return sos << "]";
}


class matrix_profile 
{
   const double MAX = std::numeric_limits<double>::max();

   double euclidean_distance (auto it_data, auto it_pattern, const std::size_t length)
   {
      double sum_sq {0.0};

      for (std::size_t i = 0; i < length; ++i)
      {
         auto diff = *(it_pattern + i) - *(it_data + i);
         sum_sq += (diff * diff);
      }
      return std::sqrt(sum_sq);
   }

   public: 
   std::vector<double> make_profile_for_1_series(const std::vector<double> & series, const std::size_t window_length)
   {
      std::vector<double> profile;
      if (series.size() < window_length)
         return profile;

      auto size = series.size();
      auto number_of_locations = size - window_length + 1;
      auto iSeries = series.begin();
      for (auto it_pattern = iSeries; it_pattern < iSeries + number_of_locations; ++it_pattern)
      {
         // dla kazdego mozliwego wzorca
         auto minimum = std::numeric_limits<double>::max();
         for (auto it_data = series.begin(); it_data < it_pattern; ++it_data)
         {
            auto distance = euclidean_distance (it_data, it_pattern, window_length);
            if (distance < minimum)
               minimum = distance;
         }
         for (auto it_data = it_pattern + 1; it_data < series.begin() + number_of_locations; ++it_data)
         {
            auto distance = euclidean_distance (it_data, it_pattern, window_length);
            if (distance < minimum)
               minimum = distance;
         }
         profile.push_back(minimum);
      }

      return profile;
   };

   std::vector<std::vector<std::vector<double>>> make_full_matrix_profile_for_multiseries(const std::vector<std::vector<double>> & multiseries, const std::size_t window_length)
   {
      std::vector<std::vector<std::vector<double>>> result;

      for (const auto & series : multiseries)
      {
         auto profile_matrix = make_full_matrix_profile_for_1_series(series, window_length);
         result.push_back(profile_matrix);
      }

      return result;
   }

   private:
   double aggregate_maximum (const std::vector<double> & values)
   {
      double result = std::numeric_limits<double>::min();
      for (const auto & val : values)
      {
         if (result < val)
            result = val;
      }
      return result;
   }

   double aggregate_vector (const std::vector<double> & values)
   {
      return aggregate_maximum(values);
   }

   private: 
   std::vector<std::vector<double>> aggregate_channels (const std::vector<std::vector<std::vector<double>>> & profiles)
   {
      std::vector<std::vector<double>> aggregation;

      auto nChannels = profiles.size();
      auto size = profiles[0].size();

      for (std::size_t row = 0; row < size; ++row)
      {
         std::vector<double> result_row;
         for (std::size_t col = 0; col < size; ++col)
         {
            std::vector<double> values;
            for (std::size_t channel = 0; channel < nChannels; ++channel)
            {
               values.push_back(profiles[channel][row][col]);
            }
            auto value = aggregate_vector (values);
            result_row.push_back(value);
         }
         // debug(row); debug(result_row);
         aggregation.push_back(result_row);
      }
      return aggregation;
   }

   private:
   double get_minimum_in_vector (const std::vector<double> & values)
   {
      double minimum = std::numeric_limits<double>::max();
      for (const auto & d : values)
      {
         if (minimum > d)
            minimum = d;
      }
      return minimum;
   }

   private:
   std::vector<double> get_minimum_for_each_row (const std::vector<std::vector<double>> & matrix)
   {
      std::vector<double> profile;
      for (const auto & row : matrix)
      {
         auto minimum = get_minimum_in_vector(row);
         profile.push_back(minimum);
      }
      return profile;
   }

   public: 
   std::vector<double> make_profile_for_multiseries (const std::vector<std::vector<double>> & time_serieses, const std::size_t window_length)
   {
      auto profiles = make_full_matrix_profile_for_multiseries(time_serieses, window_length);

      // teraz agregacja:
      auto profile_matrix = aggregate_channels(profiles);
      auto profile = get_minimum_for_each_row (profile_matrix);

      return profile;
   }

   std::vector<std::vector<double>> make_full_matrix_profile_for_1_series(const std::vector<double> & time_series, const std::size_t window_length)
   {
      std::vector<std::vector<double>> full_matrix_profile;
      if (time_series.size() < window_length)
         return full_matrix_profile;

      auto size = time_series.size();
      auto number_of_locations = size - window_length + 1;
      auto iBegin = time_series.begin();
      for (auto it_pattern = iBegin; it_pattern < iBegin + number_of_locations; ++it_pattern)
      {
         // dla kazdego mozliwego wzorca
         std::vector<double> profile;

         auto minimum = std::numeric_limits<double>::max();
         for (auto it_data = time_series.begin(); it_data < time_series.begin() + number_of_locations; ++it_data)
         {
            double distance;
            if (it_data == it_pattern)
               distance = MAX;
            else 
               distance = euclidean_distance (it_data, it_pattern, window_length);
            profile.push_back(distance);
         }
         full_matrix_profile.push_back(profile);
      }

      return full_matrix_profile;
   };
};

void exp_1 ()
{
   std::string description {"one channel simple"};
   std::cout << "=================" << std::endl;
   std::cout << __FUNCTION__ << std::endl;
   std::cout << description << std::endl;
   std::cout << "=================" << std::endl;

   std::size_t window_length {5};
   std::vector<double> series {3.4, 3.0, 3.2, 3.3, 3.1, 3.0, 3.4, 3.0, 3.2, 3.1, 7.9, 3.2, 3.5, 3.1};

   matrix_profile mp;
   auto profile = mp.make_profile_for_1_series(series, window_length);

   auto profile_matrix = mp.make_full_matrix_profile_for_1_series(series, window_length);
   debug(window_length);
   debug(series.size());
   debug(series);
   debug(profile);
   debug(profile_matrix);
}

void exp_2 ()
{
   std::string description {"one channel more complicated"};
   std::cout << "=================" << std::endl;
   std::cout << __FUNCTION__ << std::endl;
   std::cout << description << std::endl;
   std::cout << "=================" << std::endl;

   std::size_t window_length {5};
   std::vector<double> data ;

   for (double x = 0; x < 100; x += 0.1)
   {
      data.push_back(std::sin(x));
   }

   data[950] = 2;

   matrix_profile mp;
   auto profile = mp.make_profile_for_1_series(data, window_length);

   debug(window_length);
   debug(data.size());
   debug(data);
   debug(profile);
}

void exp_3 ()
{
   std::string description {"two channels"};
   std::cout << "=================" << std::endl;
   std::cout << __FUNCTION__ << std::endl;
   std::cout << description << std::endl;
   std::cout << "=================" << std::endl;

   std::size_t window_length {5};
   std::vector<double> series_1 {3.4, 3.0, 3.2, 3.3, 3.1, 3.0, 3.4, 3.0, 3.2, 3.1, 7.9, 3.2, 3.5, 3.1};

   std::vector<double> series_2 {series_1}; 
   std::rotate (series_2.begin(), series_2.begin() + (series_2.size() / 2), series_2.end());
   std::vector<std::vector<double>> multiseries { series_1, series_2 };

   matrix_profile mp;

   auto profile = mp.make_profile_for_multiseries(multiseries, window_length);
   debug(window_length);
   debug(multiseries.size());
   debug(multiseries);
   debug(profile);
}

int main()
{
   exp_1();
   // exp_2();
   exp_3();

   return 0;
}


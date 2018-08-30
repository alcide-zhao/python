# writhe the dictionary into astropytables
# from astropy.table import Table, Column
# t = Table(masked=True)
# t.add_column(Column(name='time', data=time,masked=True, unit='year'))
# t.add_column(Column(name='lon', data=lon,masked=True, unit='degree'))
# t.add_column(Column(name='lat', data=lat,masked=True, unit='lat'))
# t.add_column(Column(name='attribute', data=attribute))
# t.add_column(Column(name='RCP_jja_emissions', data=RCP_jja_emissions,unit='kg m-2 s-1'))
# Write out to file
# t.write('myfile.fits')  # also support HDF5, ASCII, etc.
 
# # Read in from file
# t = table.Table.read('myfile.fits')



# from astropy.table import Table
# t = table.Table.read('myfile.fits')
# data_rows = [(1, 2.0, 'x'),
             # (4, 5.0, 'y'),
             # (5, 8.2, 'z')]
# t = Table(rows=data_rows, names=('a', 'b', 'c'), meta={'name': 'first table'},
           # dtype=('i4', 'f8', 'S1'))
# You can also assign a unit to the columns. If any column has a unit assigned, all units would be shown as follows:
# t['b'].unit = 's'
# you can get summary information about the table as follows:
# t.info
# If you do not like the format of a particular column, you can change it
# t['b'].format = '7.3f'
# For a long table you can scroll up and down through the table one page at time:
# t.more() 
# You can also display it as an HTML-formatted table in the browser:
# t.show_in_browser()  
# or as an interactive (searchable & sortable) javascript table:
# t.show_in_browser(jsviewer=True)
# examine column names
# t.colnames
# Access the data by column or row using familiar numpy structured array syntax:
# t['a']       # Column 'a'
# You can retrieve a subset of a table by rows (using a slice) or columns (using column names), where the subset is returned as a new table:
# print(t[0:2])  
# print(t['a', 'c'])
# Replace, add, remove, and rename columns with the following:
# >>> t['b'] = ['a', 'new', 'dtype']   # Replace column b (different from in place)
# >>> t['d'] = [1, 2, 3]               # Add column d
# >>> del t['c']                       # Delete column c
# >>> t.rename_column('a', 'A')        # Rename column a to A
# Adding a new row of data to the table is as follows:
# >>> t.add_row([-8, -9, 10])
# You can create a table with support for missing values, for example by setting masked=True:
# t = Table([a, b, c], names=('a', 'b', 'c'), masked=True, dtype=('i4', 'f8', 'S1')
# You can include certain object types like Time, SkyCoord or Quantity in your table. These “mixin” columns behave like a hybrid of a regular Column and the native object type (see Mixin columns). For example:
# >>>

# >>> from astropy.time import Time
# >>> from astropy.coordinates import SkyCoord
# >>> tm = Time(['2000:002', '2002:345'])
# >>> sc = SkyCoord([10, 20], [-45, +40], unit='deg')
# >>> t = Table([tm, sc], names=['time', 'skycoord'])
# >>> t
# <Table length=2>
         # time          skycoord
                       # deg,deg
        # object          object
# --------------------- ----------
# 2000:002:00:00:00.000 10.0,-45.0
# 2002:345:00:00:00.000  20.0,40.0

/*
 * switcher.cxx
 * 
 * Copyright 2023 Mike <mike@fedora38-2.home>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * Examinr the layout/behavior of switch for handling TAG values
 * 
 */


#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
	int tag = 321;
	
	switch (tag)
	{
		case 123:
			cout << "tag 123" << endl;
			break;
		case 321:
			cout << "tag 321" << endl;
			break;
		default:
			cout << "Error: tag not known." << endl; 
			break;
	}
	return 0;
}


#pragma once


int linear_search()
{
	int hi = -1;
	for (int i = 0; i <= Sxy_size; ++i) {
		hi = (Sxy_row[i].column == Szy_element.column) ? i : hi;
	}

	if (hi >= 0) {	
		Sxy_row[hi].value += w * Sxz_elements[k].value * Szy_element.value;
	}
}


int binary_search()
{
		int lo = 0;
		int hi = Sxy_size;
								
		while (lo < hi) {
			int mid = hadd(lo, hi); // same as (lo + hi)>> 1;
			if (Sxy_row[mid].column >= Szy_element.column) {
				hi = mid;
			} else {
				lo = mid + 1;
			}
		}

		if (Sxy_row[hi].column == Szy_element.column) {	
			Sxy_row[hi].value += w * Sxz_elements[k].value * Szy_element.value;
		}
}


int binary_search_store()
{
	int lo = 0;
	int hi = Sxy_size;
	int midColumn = 0;	

	while (lo < hi) {
		int mid = hadd(lo, hi); // same as (lo + hi)>> 1;
		midColumn = Sxy_row[mid].column;
		if (midColumn >= Szy_element.column) {
			hi = mid;
		} else {
			lo = mid + 1;
		}
	}

	if (midColumn == Szy_element.column) {	
		Sxy_row[hi].value += w * Sxz_elements[k].value * Szy_element.value;
	}
}


int binary_search_unrolled()
{
	unsigned int hi = Sxy_size;
	unsigned int id = sub_sat(hi, 8u);
					
	while (Sxy_row[id].column > Szy_element.column)
		id = sub_sat(id, 16u);
	
	id = Szy_element.column < Sxy_row[id].column ? sub_sat(id, 4u) : min(id + 4, hi);
	id = Szy_element.column < Sxy_row[id].column ? sub_sat(id, 2u) : min(id + 2, hi);
	id = Szy_element.column < Sxy_row[id].column ? sub_sat(id, 1u) : min(id + 1, hi);
	id = Szy_element.column < Sxy_row[id].column ? sub_sat(id, 1u) : min(id + 0, hi);
	
	if (Sxy_row[id].column == Szy_element.column) {	
		Sxy_row[id].value += w * Sxz_elements[k].value * Szy_element.value;
		break;
	}
}


void binary_array() {
	
			// set starting position in binary search
			int initial_id = (1 << (32 - clz(Sxy_size))) >> 1; // round size up to the nearest power of 2
			int initial_step = initial_id;

			// corrections
			initial_id -= (initial_id > 0) ? 1 : 0;
			
			--Sxy_size;
			
					int id = initial_id;
					int stp = initial_step;
					int col = -1;

					while (stp) {
						stp >>= 1; // divide step by 2
						col = Sxy_row[min(id, Sxy_size)].column;								
						stp = (Szy_elem.column == col) ? 0 : stp;
						id += (Szy_elem.column < col) ? -stp : stp;
					}
								
					if (col == Szy_elem.column) {	
						Sxy_row[min(id, Sxy_size)].value += w * Sxz_elements[k].value * Szy_elem.value;
					}

}

import pytest
import hypothesis
from nexoclom.initial_state.input_classes import Geometry


@pytest.mark.initial_state
def test_geometry_search(gparam):
    geometry = Geometry(gparam)
    ids_insert = geometry.insert()
    ids_search = geometry.search()
    
    assert ids_insert == ids_search

def tile_mesh(plan, tiles):
    '''
    Construct a mesh by joining tiles according to a plan

    INPUT:
      plan: n-dimensional arrays of tile names
      tiles: dict str (tile name) -> Tile instance (that tile)

    OUTPUT:
      ntuple of (mesh, mesh functions from data)
    '''
    # Do we have all needed tiles, are they sensible?
    assert is_okay(plan, tiles)

    # Reduce 3 x 4 x 12 to
    #        3 x 4
    #        3
    #        1 tiles
    naxis = ndims(plan)
    while naxis > 0:
        plan, tiles = join_tile(plan, tiles)
        naxis -= 1

    assert len(tiles) == 1

    tile = next(iter(tiles.values()))

    

    
        
    

extends Node

func _ready():
	var voro = Voronoi.new()
	#voro.voronoi2d(5)
	#print(voro.get_faces())
	#print(voro.get_vertexes())
	#return 
	#voro.voronoi3d(25)
	#print(voro.get_faces())
	#print(voro.get_fragments())
	
	#performance_test()
	
	#var voro = Voronoi.new()
	
	voro.voronoi3d(10)
	var vertexes = voro.get_vertexes()
	var fragments = voro.get_fragments()
	var faces = voro.get_faces()

	var new_node = Spatial.new()
	new_node.name = 'FragmentNode'
	
	var UVs = PoolVector2Array()
	UVs.push_back(Vector2(0,0))
	UVs.push_back(Vector2(0,1))
	UVs.push_back(Vector2(1,1))
	UVs.push_back(Vector2(1,0))

	for counter in fragments.size():
		var face = faces[counter]
	
		var coll_shape = CollisionShape.new()
		var con_shape = ConvexPolygonShape.new()
		var body = RigidBody.new()
		body.visible = false
		body.sleeping = true
		new_node.add_child(body)
		body.add_child(coll_shape)
		body.set_owner(new_node)
		coll_shape.set_owner(new_node)
	
		var coll_points = PoolVector3Array()
		for p in face:
			coll_points.append(vertexes[p])
		con_shape.points = coll_points

		coll_shape.disabled = true
		coll_shape.name = 'CollisionShape'
		coll_shape.shape = con_shape
		var color = Color(rand_range(0.3, 0.9), rand_range(0.3, 0.9), rand_range(0.3, 0.9))
	
		var fragment = fragments[counter]
		for t in fragment:
			var tmpMesh = Mesh.new()
			var mesh = MeshInstance.new()
			var vertices = PoolVector3Array()
			var mat = SpatialMaterial.new()
			vertices.push_back(t[0])
			vertices.push_back(t[1])
			vertices.push_back(t[2])
			mat.albedo_color = color

			var st = SurfaceTool.new()
			st.begin(Mesh.PRIMITIVE_TRIANGLE_FAN)
			st.set_material(mat)
			for v in vertices.size(): 
				st.add_color(color)
				st.add_vertex(vertices[v])
				st.add_uv(UVs[v])
			st.commit(tmpMesh)
			mesh.mesh = tmpMesh
			body.add_child(mesh)
			mesh.set_owner(new_node)
	
	self.add_child(new_node)

	var packed_scene = PackedScene.new()
	packed_scene.pack(new_node)
	ResourceSaver.save("res://voronoi.tscn", packed_scene)

func performance_test():
	var voro : Voronoi
	var t

	for d in [2, 3]:
		for i in [10, 25, 50, 250, 1000]:
			voro = Voronoi.new()
			if d == 2:
				t = OS.get_ticks_msec()
				voro.voronoi2d(i)
				print('2D: ', i)
				print(OS.get_ticks_msec() - t, 'msecs')
				print('passed: ', len(voro.get_faces()) == i)
			else:
				t = OS.get_ticks_msec()
				voro.voronoi3d(i)
				print('3D: ', i)
				print(OS.get_ticks_msec() - t, 'msecs')
				print('passed: ', len(voro.get_faces()) == i and len(voro.get_fragments()) == i)
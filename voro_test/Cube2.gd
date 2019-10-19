extends KinematicBody

func _ready():
	
	var voro = Voronoi.new()
	voro.voronoi3d(100)
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
	ResourceSaver.save("res://voronoi_explosion.tscn", packed_scene)


func _physics_process(delta):
	
		
	var new_node = get_node('FragmentNode')
	new_node.transform.translated(self.transform.origin)
	get_node('CollisionShape').disabled = true
	for body in new_node.get_children():
		body.get_node('CollisionShape').disabled = false
		body.visible = true
		body.sleeping = false
			
	get_node('MeshInstance').hide()
	
func _input(event):
	if event.is_action('ui_select'):
		
		
		var new_node = get_node('FragmentNode')
		if new_node:
			new_node.free()
		var voro_cube = load('res://voronoi_explosion.tscn')
		
		voro_cube = voro_cube.instance()
		add_child(voro_cube)
		self.translate(Vector3(0, 0, 3))
		get_node('MeshInstance').show()
		get_node('CollisionShape').disabled = false;

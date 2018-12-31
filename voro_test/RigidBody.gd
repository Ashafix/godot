extends KinematicBody

var shattered = false

func _physics_process(delta):
	if shattered:
		return
	var move = move_and_collide(Vector3(0, -0.1, 0))

	if move:
		shattered = true
		print('baammm')
		var new_node = get_node('../FragmentNode')
		new_node.transform.translated(self.transform.origin)
		get_node('CollisionShape').disabled = true
		for body in new_node.get_children():
			body.get_node('CollisionShape').disabled = false
			body.visible = true
			body.sleeping = false
			
		get_node('MeshInstance').hide()

func _input(event):
	if event.is_action('ui_select') and shattered:
		print('starting again')
		shattered = false
		var new_node = get_node('../FragmentNode')
		if new_node:
			new_node.free()
		var voro_cube = load('res://voronoi.tscn')
		
		voro_cube = voro_cube.instance()
		get_parent().add_child(voro_cube)
		self.translate(Vector3(0, 12, 0))
		get_node('MeshInstance').show()
		get_node('CollisionShape').disabled = false;
		
	
	
	
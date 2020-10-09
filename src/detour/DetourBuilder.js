
class DetourBuilder {

 build( params,  tileX, tileY) {
	 console.log("NOT BUILD")
		let data = NavMeshBuilder.createNavMeshData(params);
		if (data != null) {
			data.header.x = tileX;
			data.header.y = tileY;
		}
		return data;
	}
}

export default DetourBuilder;

function Point(x, y) {
	this.x = x || 0;
	this.y = y || 0;

	this.__defineGetter__('length', function() {
		return Math.sqrt(this.x * this.x + this.y * this.y);
	});

	this.clone = function() {
		return new Point(this.x, this.y);
	};

	this.add = function(point) {
		return new Point(this.x + point.x, this.y + point.y);
	};

	this.subtract = function(point) {
		return new Point(point.x - this.x, point.y - this.y);
	};

	this.offset = function(value) {
		this.x += value;
		this.y += value;
	};

	this.distance = function(x2, y2) {
		var dx = x2 - this.x;
		var dy = y2 - this.y;
		return Math.sqrt(dx * dx + dy * dy);
	};

	this.equals = function(point) {
		return this.x === point.x && this.y === point.y;
	};

	this.normalize = function(thickness) {
		var ratio = thickness / this.length;
		this.x += this.x * ratio;
		this.y += this.y * ratio;
	};

	this.interpolate = function(p1, p2, f) {
		var point = new Point();
		point.x = p1.x + f * (p2.x - p1.x);
		point.y = p1.y + f * (p2.y - p1.y);
		return point;
	};

	this.polar = function(len, angle) {
		this.x = Math.cos(angle) * len;
		this.y = Math.sin(angle) * len;
	}
}
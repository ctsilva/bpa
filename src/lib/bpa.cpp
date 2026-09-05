#include "bpa.h"

#include <algorithm>
#include <cstdint>
#include <deque>
#include <limits>
#include <numbers>
#include <numeric>
#include <optional>
#include <vector>

using namespace glm;

namespace bpa {
	namespace {
		struct MeshEdge;

		struct MeshPoint {
			dvec3 pos;
			dvec3 normal;
			std::uint32_t index; // into the input vector
			bool used = false;
			std::vector<MeshEdge*> edges;
		};

		enum class EdgeStatus { active, inner, boundary };

		struct MeshEdge {
			MeshPoint* a;
			MeshPoint* b;
			MeshPoint* opposite;
			dvec3 center;
			MeshEdge* prev;
			MeshEdge* next;
			EdgeStatus status = EdgeStatus::active;
		};

		struct MeshFace : std::array<MeshPoint*, 3> {
			auto normal() const { return normalize(cross((*this)[0]->pos - (*this)[1]->pos, (*this)[0]->pos - (*this)[2]->pos)); }
		};

		using Cell = std::vector<MeshPoint>;

		// Uniform grid of cells of side 2 * radius: every point a ball touching a point in a cell
		// can touch lies in that cell or one of its 26 neighbours.
		struct Grid {
			Grid(const std::vector<Point>& points, double radius) : cellSize(radius * 2) {
				lower = points.front().pos;
				upper = points.front().pos;
				for (const auto& p : points) {
					lower = min(lower, p.pos);
					upper = max(upper, p.pos);
				}
				dims = max(ivec3{ceil((upper - lower) / cellSize)}, ivec3{1});
				cells.resize(std::size_t(dims.x) * dims.y * dims.z);
				for (std::size_t i = 0; i < points.size(); i++)
					cell(cellIndex(points[i].pos)).push_back({points[i].pos, points[i].normal, static_cast<std::uint32_t>(i)});
			}

			auto cellIndex(dvec3 point) const -> ivec3 { return clamp(ivec3{(point - lower) / cellSize}, ivec3{}, dims - 1); }

			auto cell(ivec3 index) -> Cell& { return cells[std::size_t(index.z) * dims.x * dims.y + std::size_t(index.y) * dims.x + index.x]; }

			// The points within 2 * radius of `point`, except those in `ignore`.
			auto sphericalNeighborhood(dvec3 point, std::initializer_list<const MeshPoint*> ignore) -> std::vector<MeshPoint*> {
				std::vector<MeshPoint*> result;
				const auto centerIndex = cellIndex(point);
				result.reserve(cell(centerIndex).size() * 27);
				for (auto xOff : {-1, 0, 1})
					for (auto yOff : {-1, 0, 1})
						for (auto zOff : {-1, 0, 1}) {
							const auto index = centerIndex + ivec3{xOff, yOff, zOff};
							if (any(lessThan(index, ivec3{})) || any(greaterThanEqual(index, dims)))
								continue;
							for (auto& p : cell(index))
								if (length2(p.pos - point) < cellSize * cellSize && std::find(begin(ignore), end(ignore), &p) == end(ignore))
									result.push_back(&p);
						}
				return result;
			}

			dvec3 lower;
			dvec3 upper;
			double cellSize;
			ivec3 dims;
			std::vector<Cell> cells;
		};

		// Centre of the ball of the given radius through the three points of `f`, on the side of
		// its normal; nothing if the circumradius is larger than the ball.
		auto computeBallCenter(MeshFace f, double radius) -> std::optional<dvec3> {
			const dvec3 ac = f[2]->pos - f[0]->pos;
			const dvec3 ab = f[1]->pos - f[0]->pos;
			const dvec3 abXac = cross(ab, ac);
			const auto nn = dot(abXac, abXac);
			if (nn <= 1e-20 * dot(ab, ab) * dot(ac, ac)) // degenerate triangle
				return {};
			const dvec3 toCircumCircleCenter = (cross(abXac, ab) * dot(ac, ac) + cross(ac, abXac) * dot(ab, ab)) / (2 * nn);
			const dvec3 circumCircleCenter = f[0]->pos + toCircumCircleCenter;

			const auto heightSquared = radius * radius - dot(toCircumCircleCenter, toCircumCircleCenter);
			if (heightSquared < 0)
				return {};
			return circumCircleCenter + f.normal() * std::sqrt(heightSquared);
		}

		// Points within a relative 1e-9 of the sphere count as outside, so that the three points
		// the ball touches, and points exactly cospherical with them, do not fail the test.
		constexpr auto emptinessTolerance = 1e-9;

		auto ballIsEmpty(dvec3 ballCenter, const std::vector<MeshPoint*>& points, double radius) -> bool {
			const auto r2 = radius * (1 - emptinessTolerance) * radius * (1 - emptinessTolerance);
			return !std::any_of(begin(points), end(points), [&](MeshPoint* p) { return length2(p->pos - ballCenter) < r2; });
		}

		struct SeedResult {
			MeshFace f;
			dvec3 ballCenter;
		};

		// Only the nearest this many neighbours of a seed candidate are paired: a valid seed's
		// other two vertices are almost always among the closest points, and for a candidate
		// under an already reconstructed sheet every pair fails, so the pair loop must be cheap.
		constexpr std::size_t seedNeighbors = 100;

		// Seed search (section 4.3 of the paper). Cells are visited from a cursor that persists
		// between calls; a cell holding a used point is skipped (the paper's heuristic against
		// spawning small components next to the surface, fig. 4c), else one candidate is tried,
		// the point projecting furthest along the cell's average normal, paired with its nearest
		// neighbours: the first triangle with an empty ball on the side of that normal is the
		// seed. A candidate that fails cannot succeed later (points only become used), so the
		// cursor never moves back; after a seed it stays, and the cell is skipped next time.
		auto findSeedTriangle(Grid& grid, double radius, std::size_t& cursor) -> std::optional<SeedResult> {
			for (; cursor < grid.cells.size(); cursor++) {
				auto& cell = grid.cells[cursor];
				if (cell.empty() || std::any_of(begin(cell), end(cell), [](const MeshPoint& p) { return p.used; }))
					continue;
				const auto avgNormal = normalize(std::accumulate(begin(cell), end(cell), dvec3{}, [](dvec3 acc, const MeshPoint& p) { return acc + p.normal; }));
				const auto centroid = std::accumulate(begin(cell), end(cell), dvec3{}, [](dvec3 acc, const MeshPoint& p) { return acc + p.pos; }) / double(cell.size());
				auto& p1 = *std::max_element(begin(cell), end(cell), [&](const MeshPoint& a, const MeshPoint& b) {
					return dot(a.pos - centroid, avgNormal) < dot(b.pos - centroid, avgNormal);
				});

				auto neighborhood = grid.sphericalNeighborhood(p1.pos, {&p1});
				std::sort(begin(neighborhood), end(neighborhood), [&](MeshPoint* a, MeshPoint* b) { return length(a->pos - p1.pos) < length(b->pos - p1.pos); });
				const auto nPairs = std::min(neighborhood.size(), seedNeighbors);
				for (std::size_t i2 = 0; i2 < nPairs; i2++) {
					auto* p2 = neighborhood[i2];
					if (p2->used)
						continue;
					for (std::size_t i3 = 0; i3 < nPairs; i3++) {
						auto* p3 = neighborhood[i3];
						if (p2 == p3 || p3->used)
							continue;
						MeshFace f{{&p1, p2, p3}};
						if (dot(f.normal(), avgNormal) < 0)
							continue;
						const auto ballCenter = computeBallCenter(f, radius);
						if (ballCenter && ballIsEmpty(ballCenter.value(), neighborhood, radius)) {
							p1.used = true;
							p2->used = true;
							p3->used = true;
							return SeedResult{f, ballCenter.value()};
						}
					}
				}
			}
			return {};
		}

		auto getActiveEdge(std::vector<MeshEdge*>& front) -> std::optional<MeshEdge*> {
			while (!front.empty()) {
				auto* e = front.back();
				if (e->status == EdgeStatus::active)
					return e;
				front.pop_back(); // non-active edges are dropped lazily
			}
			return {};
		}

		struct PivotResult {
			MeshPoint* p;
			dvec3 center;
		};

		// Roll the ball around edge `e`, starting from its stored centre, and return the point it
		// touches first together with the centre at that moment; nothing if it touches no point
		// or the ball there is not empty.
		auto ballPivot(const MeshEdge* e, Grid& grid, double radius) -> std::optional<PivotResult> {
			const auto m = (e->a->pos + e->b->pos) / 2.0;
			const auto oldCenterVec = normalize(e->center - m);
			auto neighborhood = grid.sphericalNeighborhood(m, {e->a, e->b, e->opposite});

			auto smallestAngle = std::numeric_limits<double>::max();
			MeshPoint* pointWithSmallestAngle = nullptr;
			dvec3 centerOfSmallest{};
			for (const auto& p : neighborhood) {
				const auto newFaceNormal = MeshFace{{e->b, e->a, p}}.normal();

				// not in the paper: the new triangle must agree with the normal of the point
				if (dot(newFaceNormal, p->normal) < 0)
					continue;

				const auto c = computeBallCenter(MeshFace{{e->b, e->a, p}}, radius);
				if (!c)
					continue;

				// not in the paper: the ball centre must be above the new triangle
				const auto newCenterVec = normalize(c.value() - m);
				if (dot(newCenterVec, newFaceNormal) < 0)
					continue;

				// not in the paper: points joined to an edge end by an inner edge are skipped
				const auto innerEdgeExists = std::any_of(begin(p->edges), end(p->edges), [&](const MeshEdge* ee) {
					const auto* otherPoint = ee->a == p ? ee->b : ee->a;
					return ee->status == EdgeStatus::inner && (otherPoint == e->a || otherPoint == e->b);
				});
				if (innerEdgeExists)
					continue;

				auto angle = std::acos(std::clamp(dot(oldCenterVec, newCenterVec), -1.0, 1.0));
				if (dot(cross(newCenterVec, oldCenterVec), e->a->pos - e->b->pos) < 0)
					angle += std::numbers::pi;
				if (angle < smallestAngle) {
					smallestAngle = angle;
					pointWithSmallestAngle = p;
					centerOfSmallest = c.value();
				}
			}

			if (pointWithSmallestAngle && ballIsEmpty(centerOfSmallest, neighborhood, radius))
				return PivotResult{pointWithSmallestAngle, centerOfSmallest};
			return {};
		}

		auto notUsed(const MeshPoint* p) -> bool { return !p->used; }

		auto onFront(const MeshPoint* p) -> bool {
			return std::any_of(begin(p->edges), end(p->edges), [&](const MeshEdge* e) { return e->status == EdgeStatus::active; });
		}

		// Mark the edge inner; getActiveEdge() drops it from the front later.
		void remove(MeshEdge* edge) { edge->status = EdgeStatus::inner; }

		void outputTriangle(MeshFace f, std::vector<Face>& triangles) { triangles.push_back({f[0]->index, f[1]->index, f[2]->index}); }

		auto join(MeshEdge* e_ij, MeshPoint* o_k, dvec3 o_k_ballCenter, std::vector<MeshEdge*>& front, std::deque<MeshEdge>& edges)
			-> std::pair<MeshEdge*, MeshEdge*> {
			auto& e_ik = edges.emplace_back(MeshEdge{e_ij->a, o_k, e_ij->b, o_k_ballCenter});
			auto& e_kj = edges.emplace_back(MeshEdge{o_k, e_ij->b, e_ij->a, o_k_ballCenter});

			e_ik.next = &e_kj;
			e_ik.prev = e_ij->prev;
			e_ij->prev->next = &e_ik;
			e_ij->a->edges.push_back(&e_ik);

			e_kj.prev = &e_ik;
			e_kj.next = e_ij->next;
			e_ij->next->prev = &e_kj;
			e_ij->b->edges.push_back(&e_kj);

			o_k->used = true;
			o_k->edges.push_back(&e_ik);
			o_k->edges.push_back(&e_kj);

			front.push_back(&e_ik);
			front.push_back(&e_kj);
			remove(e_ij);

			return {&e_ik, &e_kj};
		}

		void glue(MeshEdge* a, MeshEdge* b) {
			// case 1: a and b form a loop of two
			if (a->next == b && a->prev == b && b->next == a && b->prev == a) {
				remove(a);
				remove(b);
				return;
			}
			// case 2: a and b are adjacent
			if (a->next == b && b->prev == a) {
				a->prev->next = b->next;
				b->next->prev = a->prev;
				remove(a);
				remove(b);
				return;
			}
			if (a->prev == b && b->next == a) {
				a->next->prev = b->prev;
				b->prev->next = a->next;
				remove(a);
				remove(b);
				return;
			}
			// cases 3 and 4: a and b are apart, in the same loop or in two
			a->prev->next = b->next;
			b->next->prev = a->prev;
			a->next->prev = b->prev;
			b->prev->next = a->next;
			remove(a);
			remove(b);
		}

		auto findReverseEdgeOnFront(MeshEdge* edge) -> MeshEdge* {
			for (auto& e : edge->a->edges)
				if (e->a == edge->b)
					return e;
			return nullptr;
		}
	} // namespace

	auto reconstruct(const std::vector<Point>& points, double radius) -> std::vector<Face> {
		if (points.empty())
			return {};
		Grid grid(points, radius);

		std::vector<Face> triangles;
		std::deque<MeshEdge> edges; // stable addresses
		std::vector<MeshEdge*> front;
		std::size_t seedCursor = 0;

		// Fig. 5 of the paper: pivot until the front is exhausted, seed again among the points
		// still unused, until no seed is left.
		while (const auto seedResult = findSeedTriangle(grid, radius, seedCursor)) {
			auto [seed, ballCenter] = seedResult.value();
			outputTriangle(seed, triangles);
			auto& e0 = edges.emplace_back(MeshEdge{seed[0], seed[1], seed[2], ballCenter});
			auto& e1 = edges.emplace_back(MeshEdge{seed[1], seed[2], seed[0], ballCenter});
			auto& e2 = edges.emplace_back(MeshEdge{seed[2], seed[0], seed[1], ballCenter});
			e0.prev = e1.next = &e2;
			e0.next = e2.prev = &e1;
			e1.prev = e2.next = &e0;
			seed[0]->edges = {&e0, &e2};
			seed[1]->edges = {&e0, &e1};
			seed[2]->edges = {&e1, &e2};
			front.insert(end(front), {&e0, &e1, &e2});

			while (auto e_ij = getActiveEdge(front)) {
				const auto o_k = ballPivot(e_ij.value(), grid, radius);
				if (o_k && (notUsed(o_k->p) || onFront(o_k->p))) {
					outputTriangle({{e_ij.value()->a, o_k->p, e_ij.value()->b}}, triangles);
					auto [e_ik, e_kj] = join(e_ij.value(), o_k->p, o_k->center, front, edges);
					if (auto* e_ki = findReverseEdgeOnFront(e_ik))
						glue(e_ik, e_ki);
					if (auto* e_jk = findReverseEdgeOnFront(e_kj))
						glue(e_kj, e_jk);
				} else {
					e_ij.value()->status = EdgeStatus::boundary;
				}
			}
		}

		return triangles;
	}
} // namespace bpa

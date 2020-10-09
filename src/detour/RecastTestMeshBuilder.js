/*
Recast4J Copyright (c) 2015 Piotr Piastucki piotr@jtilia.org

This software is provided 'as-is', without any express or implied
warranty.  In no event will the authors be held liable for any damages
arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:
1. The origin of this software must not be misrepresented; you must not
 claim that you wrote the original software. If you use this software
 in a product, an acknowledgment in the product documentation would be
 appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
 misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
import NavMesh from "./NavMesh.js"
import ObjImporter from "../recast/ObjImporter.js"
import RecastConstants from "../recast/RecastConstants.js"
import RecastConfig from "../recast/RecastConfig.js"
import SampleAreaModifications from "../detour/SampleAreaModifications.js"
import RecastBuilderConfig from "../recast/RecastBuilderConfig.js"
import RecastBuilder from "../recast/RecastBuilder.js"
import NavMeshDataCreateParams from "./NavMeshDataCreateParams.js"
import NavMeshBuilder from "./NavMeshBuilder.js"

class RecastTestMeshBuilder {

    meshData;
    static m_cellSize = 0.3;
    static m_cellHeight = 0.2;
    static m_agentHeight = 2.0;
    static m_agentRadius = 0.6;
    static m_agentMaxClimb = 0.9;
    static m_agentMaxSlope = 45.0;
    static m_regionMinSize = 8;
    static m_regionMergeSize = 20;
    static m_edgeMaxLen = 12.0;
    static m_edgeMaxError = 1.3;
    static m_vertsPerPoly = 6;
    static m_detailSampleDist = 6.0;
    static m_detailSampleMaxError = 1.0;

    static fromFile = function(objFileContents) {
        return new RecastTestMeshBuilder(new ObjImporter().load(objFileContents), RecastConstants.WATERSHED,
            this.m_cellSize, this.m_cellHeight, this.m_agentHeight,this.m_agentRadius, this.m_agentMaxClimb, this.m_agentMaxSlope,
            this.m_regionMinSize, this.m_regionMergeSize, this.m_edgeMaxLen, this.m_edgeMaxError, this.m_vertsPerPoly, this.m_detailSampleDist,
            this.m_detailSampleMaxError);
    }

    constructor(m_geom, m_partitionType, m_cellSize, m_cellHeight, m_agentHeight, m_agentRadius, m_agentMaxClimb, m_agentMaxSlope, m_regionMinSize,
        m_regionMergeSize, m_edgeMaxLen, m_edgeMaxError,  m_vertsPerPoly,
        m_detailSampleDist, m_detailSampleMaxError) {
         let cfg = new RecastConfig(m_partitionType, m_cellSize, m_cellHeight, m_agentHeight, m_agentRadius,
            m_agentMaxClimb, m_agentMaxSlope, m_regionMinSize, m_regionMergeSize, m_edgeMaxLen, m_edgeMaxError,
            m_vertsPerPoly, m_detailSampleDist, m_detailSampleMaxError, 0, SampleAreaModifications.SAMPLE_AREAMOD_GROUND);
        let bcfg = new RecastBuilderConfig(cfg, m_geom.getMeshBoundsMin(), m_geom.getMeshBoundsMax());
        let rcBuilder = new RecastBuilder();
        let rcResult = rcBuilder.build(m_geom, bcfg);
        let m_pmesh = rcResult.getMesh();
        for (let i = 0; i < m_pmesh.npolys; ++i) {
            m_pmesh.flags[i] = 1;
        }
        let m_dmesh = rcResult.getMeshDetail();
        let params = new NavMeshDataCreateParams();
        params.verts = m_pmesh.verts;
        params.vertCount = m_pmesh.nverts;
        params.polys = m_pmesh.polys;
        params.polyAreas = m_pmesh.areas;
        params.polyFlags = m_pmesh.flags;
        params.polyCount = m_pmesh.npolys;
        params.nvp = m_pmesh.nvp;
        params.detailMeshes = m_dmesh.meshes;
        params.detailVerts = m_dmesh.verts;
        params.detailVertsCount = m_dmesh.nverts;
        params.detailTris = m_dmesh.tris;
        params.detailTriCount = m_dmesh.ntris;
        params.walkableHeight = m_agentHeight;
        params.walkableRadius = m_agentRadius;
        params.walkableClimb = m_agentMaxClimb;
        params.bmin = m_pmesh.bmin;
        params.bmax = m_pmesh.bmax;
        params.cs = m_cellSize;
        params.ch = m_cellHeight;
        params.buildBvTree = true;

        params.offMeshConVerts = new Array(6);
        params.offMeshConVerts[0] = 0.1;
        params.offMeshConVerts[1] = 0.2;
        params.offMeshConVerts[2] = 0.3;
        params.offMeshConVerts[3] = 0.4;
        params.offMeshConVerts[4] = 0.5;
        params.offMeshConVerts[5] = 0.6;
        params.offMeshConRad = new Array(1);
        params.offMeshConRad[0] = 0.1;
        params.offMeshConDir = new Array(1);
        params.offMeshConDir[0] = 1;
        params.offMeshConAreas = new Array(1);
        params.offMeshConAreas[0] = 2;
        params.offMeshConFlags = new Array(1);
        params.offMeshConFlags[0] = 12;
        params.offMeshConUserID = new Array(1);
        params.offMeshConUserID[0] = 0x4567;
        params.offMeshConCount = 1;
        this.meshData = NavMeshBuilder.createNavMeshData(params);
    }

    getMeshData() {
        return this.meshData;
    }
}

export default RecastTestMeshBuilder